#!/usr/bin/env python3
"""20GA0316-centered sol HGT analysis. No HGT conclusion is hard-coded.
Usage: python 03_sol_HGT.py discover|phylogeny --root /scratch/al98750/Roux
"""
import argparse
import csv
import json
import math
import random
import re
import shutil
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from Bio import AlignIO, Phylo, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

PRIMARY = 'GCF_020740305.1'
SOURCE = 'GCA_020740305.1'
SOL_TAGS = dict(zip(['solA','solB','solC','solD','solE','solF','solG','solH','solI','solM','solJ','solK','solL'],
                    [f'LLR01_{i}' for i in range(11590,11529,-5)]))
SOL_TREE_GENES = ['solB','solC','solD','solE','solI','solK','solL']

def canonical_contig(name):
    return name.split('.')[0].removeprefix('NZ_')
EXTRA = {'GCA_000365285.1':'Dickeya solani MK10',
         'GCA_000295955':'Pantoea sp. A4',
         'GCA_900068895':'Erwinia sp. ErVv1'}


def run(cmd, stdout=None, cwd=None):
    print('+', ' '.join(map(str,cmd)), flush=True)
    subprocess.run([str(x) for x in cmd], check=True, stdout=stdout, cwd=cwd)


def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def read_tsv(path):
    with path.open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))


def complete(work, step, fn, checks=()):
    flag=work/'.state'/f'{step}.done'
    if flag.exists() and all(Path(p).is_file() and Path(p).stat().st_size>0 for p in checks):
        print(f'SKIP {step} (checkpoint)', flush=True); return
    print(f'RUN {step}',flush=True)
    fn()
    for p in checks:
        if not Path(p).is_file() or Path(p).stat().st_size == 0:
            raise RuntimeError(f'{step}: expected nonempty output missing: {p}')
    flag.write_text('success\n')


def setup(root):
    work=root/'07_sol_HGT'
    for d in ['.state','01_genomes','02_reference','03_proteomes','04_blast','05_results',
              '06_neighborhoods','07_phylogeny','08_composition','09_flank_anchors','tmp']:
        (work/d).mkdir(parents=True,exist_ok=True)
    return work


def download(acc, dest, includes='genome'):
    dest.mkdir(parents=True,exist_ok=True)
    z=dest/'package.zip'
    if not z.exists() or z.stat().st_size==0:
        run(['datasets','download','genome','accession',acc,'--include',includes,'--filename',z])
    with zipfile.ZipFile(z) as archive:
        archive.extractall(dest/'package')
    return dest/'package'/'ncbi_dataset'/'data'


def pick_genome(directory):
    hits=sorted(directory.rglob('*.fna'))
    hits=[p for p in hits if 'genomic' in p.name] or hits
    if not hits: raise RuntimeError(f'No genomic FASTA found in {directory}')
    return hits[0]


def discover_genomes(root,work):
    rows=[]; selected=set()
    for f in sorted((root/'02_genomes').glob('*.fna')):
        acc=f.stem
        if acc not in selected:
            dest=work/'01_genomes'/f'{acc}.fna'
            if not dest.exists(): dest.symlink_to(f.resolve())
            rows.append({'accession':acc,'organism':'Rouxiella badensis','source':'existing_Roux'})
            selected.add(acc)
    if PRIMARY not in selected: raise RuntimeError(f'Missing primary genome: {root}/02_genomes/{PRIMARY}.fna')
    # Retrieve Rouxiella genus; keep up to three additional genomes per non-badensis species.
    tsv=work/'tmp'/'rouxiella_ncbi.tsv'
    if not tsv.exists() or tsv.stat().st_size==0:
        p1=subprocess.Popen(['datasets','summary','genome','taxon','Rouxiella'],stdout=subprocess.PIPE)
        with tsv.open('w') as o:
            p2=subprocess.run(['dataformat','tsv','genome','--fields','accession,organism-name'],stdin=p1.stdout,stdout=o)
        p1.stdout.close(); p1.wait()
        if p1.returncode or p2.returncode: raise RuntimeError('NCBI Rouxiella genus discovery failed')
    by_species=defaultdict(list)
    for r in read_tsv(tsv):
        acc=r.get('Assembly Accession') or r.get('accession') or next((v for v in r.values() if re.match(r'^GC[AF]_',v)),None)
        name=r.get('Organism Name') or r.get('organism-name') or next((v for v in r.values() if 'Rouxiella ' in v),'')
        if not acc or not name.startswith('Rouxiella '): continue
        species=' '.join(name.split()[:2]); by_species[species].append((acc,name))
    for species, candidates in sorted(by_species.items()):
        if species=='Rouxiella badensis': continue
        # De-duplicate GCA/GCF paired assemblies by numeric accession; prefer RefSeq.
        uniq={}
        for acc,name in candidates:
            key=acc.split('_',1)[-1].split('.')[0]
            if key not in uniq or acc.startswith('GCF_'): uniq[key]=(acc,name)
        for acc,name in sorted(uniq.values())[:3]:
            if acc not in selected:
                selected.add(acc); rows.append({'accession':acc,'organism':name,'source':'NCBI_Rouxiella'})
    for acc,name in EXTRA.items():
        if acc not in selected:
            selected.add(acc); rows.append({'accession':acc,'organism':name,'source':'NCBI_outgroup_sol_positive'})
    # Download one assembly at a time to preserve retryability.
    for r in rows:
        acc=r['accession']; dest=work/'01_genomes'/f'{acc}.fna'
        if dest.exists() and dest.stat().st_size>0: continue
        pack=download(acc,work/'tmp'/'genome_packages'/acc,'genome,gbff')
        shutil.copy2(pick_genome(pack),dest)
        gbffs=list(pack.rglob('*.gbff'))
        if gbffs: shutil.copy2(gbffs[0],work/'01_genomes'/f'{acc}.gbff')
    write_tsv(work/'05_results'/'genome_manifest.tsv',rows,['accession','organism','source'])
    if len({r['organism'].split()[1] for r in rows if r['organism'].startswith('Rouxiella ')})<2:
        raise RuntimeError('Only one Rouxiella species retrieved; need genus-level comparison. Check NCBI output.')
    print(f'Genomes ready: {len(rows)}')


def extract_reference(work):
    pack=download(SOURCE,work/'02_reference'/'source','genome,gbff,protein')
    gbffs=list(pack.rglob('*.gbff'))
    if not gbffs: raise RuntimeError('GenBank source GBFF unavailable: cannot validate original LLR01 sol locus tags')
    found={}; reference_contig=None; sequence=None
    for rec in SeqIO.parse(gbffs[0],'genbank'):
        for feature in rec.features:
            if feature.type!='CDS': continue
            tag=feature.qualifiers.get('locus_tag',[''])[0]
            for gene,expected in SOL_TAGS.items():
                if tag==expected:
                    if gene in found: raise RuntimeError(f'Duplicate locus tag {tag}')
                    aa=feature.qualifiers.get('translation',[''])[0]
                    if not aa: aa=str(feature.extract(rec.seq).translate(table=11,to_stop=True))
                    found[gene]=(rec.id,int(feature.location.start),int(feature.location.end),aa)
                    if reference_contig and reference_contig!=rec.id:
                        raise RuntimeError('sol genes are not all on one contig')
                    reference_contig=rec.id; sequence=rec.seq
    missing=sorted(set(SOL_TAGS)-set(found))
    if missing: raise RuntimeError(f'Missing source sol locus tags: {missing}; inspect GBFF and do not infer absent genes')
    if not canonical_contig(reference_contig).startswith('JAJGAU010000005'):
        raise RuntimeError(f'Unexpected sol contig: {reference_contig}; expected JAJGAU010000005')
    with (work/'02_reference'/'sol_queries.faa').open('w') as o:
        for gene in SOL_TAGS:
            o.write(f'>{gene} {SOL_TAGS[gene]} 20GA0316\n{found[gene][3]}\n')
    start=min(x[1] for x in found.values()); end=max(x[2] for x in found.values())
    data={'contig':reference_contig,'start_0based':start,'end_exclusive':end,
          'length':end-start,'genes':{k:{'tag':SOL_TAGS[k],'start_0based':v[1],'end_exclusive':v[2]} for k,v in found.items()}}
    (work/'02_reference'/'sol_reference.json').write_text(json.dumps(data,indent=2)+'\n')
    (work/'02_reference'/'sol_reference.fna').write_text(f'>20GA0316_sol_region\n{sequence[start:end]}\n')
    print(f'Validated {len(found)} sol genes; {reference_contig}:{start+1}-{end} ({end-start} bp)')


def predict(work):
    manifest=read_tsv(work/'05_results'/'genome_manifest.tsv')
    for r in manifest:
        acc=r['accession']; fna=work/'01_genomes'/f'{acc}.fna'
        faa=work/'03_proteomes'/f'{acc}.faa'; gff=work/'03_proteomes'/f'{acc}.gff'
        flag=work/'.state'/f'predict_{acc}.done'
        if flag.exists() and faa.exists() and faa.stat().st_size and gff.exists() and gff.stat().st_size: continue
        run(['prodigal','-i',fna,'-a',faa,'-o',gff,'-f','gff','-p','single','-q'])
        if not faa.exists() or faa.stat().st_size==0: raise RuntimeError(f'No proteins: {acc}')
        flag.write_text('success\n')


def proteins(work,acc):
    """Return protein ID -> (contig, start0, end, strand, sequence)."""
    out={}
    fna_ids=set(x.id for x in SeqIO.parse(work/'01_genomes'/f'{acc}.fna','fasta'))
    for rec in SeqIO.parse(work/'03_proteomes'/f'{acc}.faa','fasta'):
        parts=rec.description.split(' # ')
        if len(parts)<4: raise RuntimeError(f'Unexpected Prodigal header: {rec.description}')
        contig=rec.id.rsplit('_',1)[0]
        if contig not in fna_ids: raise RuntimeError(f'Prodigal contig mismatch {contig} in {acc}')
        out[rec.id]=(contig,int(parts[1])-1,int(parts[2]),int(parts[3]),str(rec.seq))
    return out


def blast(work):
    query=work/'02_reference'/'sol_queries.faa'
    for r in read_tsv(work/'05_results'/'genome_manifest.tsv'):
        acc=r['accession']; faa=work/'03_proteomes'/f'{acc}.faa'
        db=work/'04_blast'/f'{acc}_db'; hits=work/'04_blast'/f'{acc}.tsv'
        flag=work/'.state'/f'blast_{acc}.done'
        if flag.exists() and hits.exists(): continue # Empty hits valid after successful BLAST.
        run(['makeblastdb','-in',faa,'-dbtype','prot','-out',db],stdout=subprocess.DEVNULL)
        with hits.open('w') as out:
            run(['blastp','-query',query,'-db',db,'-evalue','1e-5','-max_target_seqs','100',
                 '-num_threads','4','-outfmt','6 qseqid sseqid pident length qlen slen qcovs evalue bitscore'],stdout=out)
        flag.write_text('success\n')


def get_hits(work,acc,protein_map):
    best={}
    with (work/'04_blast'/f'{acc}.tsv').open() as f:
        for line in f:
            a=line.split();
            if len(a)!=9: continue
            q,s=a[:2]
            if q not in SOL_TAGS or s not in protein_map: continue
            ident=float(a[2]); cov=float(a[6]); ev=float(a[7]); bits=float(a[8])
            if ident<30 or cov<40 or ev>1e-5: continue
            key=(q,s)
            if key not in best or bits>best[key]['bitscore']:
                contig,start,end,strand,_=protein_map[s]
                best[key]=dict(gene=q,protein=s,contig=contig,start=start,end=end,strand=strand,
                               identity=ident,coverage=cov,evalue=ev,bitscore=bits)
    return list(best.values())


def select_locus(hits,window=85000):
    """Find a compact, multi-gene candidate; prevent one protein counting as multiple sol genes."""
    by_contig=defaultdict(list)
    for h in hits: by_contig[h['contig']].append(h)
    candidates=[]
    for contig,arr in by_contig.items():
        arr.sort(key=lambda x:x['start'])
        for anchor in arr:
            within=[h for h in arr if h['start']>=anchor['start'] and h['start']-anchor['start']<=window]
            chosen={}; used=set()
            for h in sorted(within,key=lambda x:x['bitscore'],reverse=True):
                if h['gene'] not in chosen and h['protein'] not in used:
                    chosen[h['gene']]=h; used.add(h['protein'])
            if not chosen: continue
            vals=list(chosen.values()); span=max(h['end'] for h in vals)-min(h['start'] for h in vals)
            candidates.append((len(chosen),sum(h['bitscore'] for h in vals),-span,contig,chosen))
    if not candidates: return None,{}
    top=max(candidates,key=lambda x:x[:3])
    return top[3],top[4]


def summarize(work):
    manifest=read_tsv(work/'05_results'/'genome_manifest.tsv')
    locus_rows=[]; gene_rows=[]; matrix=[]; mapping={}
    for r in manifest:
        acc=r['accession']; pmap=proteins(work,acc); hits=get_hits(work,acc,pmap)
        contig,chosen=select_locus(hits)
        n=len(chosen); span=max((h['end'] for h in chosen.values()),default=0)-min((h['start'] for h in chosen.values()),default=0)
        if n>=9 and span<=85000: status='candidate_intact_cluster'
        elif n>=4: status='partial_or_ambiguous'
        else: status='no_coherent_cluster_detected'
        seqs=SeqIO.to_dict(SeqIO.parse(work/'01_genomes'/f'{acc}.fna','fasta'))
        left=min((h['start'] for h in chosen.values()),default=0)
        right=max((h['end'] for h in chosen.values()),default=0)
        edge=bool(contig and (left<10000 or len(seqs[contig])-right<10000))
        locus_rows.append(dict(accession=acc,organism=r['organism'],status=status,genes_detected=n,
                               contig=contig or 'NA',start_1based=left+1 if contig else 'NA',
                               end=right if contig else 'NA',span_bp=span if contig else 'NA',
                               near_contig_edge=edge,genome_contigs=len(seqs)))
        matrix.append(dict(accession=acc,organism=r['organism'],status=status,**{g:int(g in chosen) for g in SOL_TAGS}))
        for gene,h in chosen.items(): gene_rows.append(dict(accession=acc,**h))
        mapping[acc]=chosen
        if acc==PRIMARY and n<10:
            raise RuntimeError(f'Primary 20GA0316 yielded only {n} sol proteins; check reference and Prodigal before proceeding')
    write_tsv(work/'05_results'/'sol_loci.tsv',locus_rows,
              ['accession','organism','status','genes_detected','contig','start_1based','end','span_bp','near_contig_edge','genome_contigs'])
    write_tsv(work/'05_results'/'sol_gene_matrix.tsv',matrix,['accession','organism','status']+list(SOL_TAGS))
    write_tsv(work/'05_results'/'sol_gene_hits.tsv',gene_rows,
              ['accession','gene','protein','contig','start','end','strand','identity','coverage','evalue','bitscore'])
    return mapping


def locus_mapping(work):
    d=defaultdict(dict)
    for r in read_tsv(work/'05_results'/'sol_gene_hits.tsv'):
        for k in ['start','end','strand']:r[k]=int(r[k])
        for k in ['identity','coverage','evalue','bitscore']:r[k]=float(r[k])
        d[r['accession']][r['gene']]=r
    return d


def neighborhoods(work):
    loci={r['accession']:r for r in read_tsv(work/'05_results'/'sol_loci.tsv')}
    mapping=locus_mapping(work)
    for acc,row in loci.items():
        if row['status']!='candidate_intact_cluster': continue
        hits=mapping[acc]; contig=row['contig']
        seq=SeqIO.to_dict(SeqIO.parse(work/'01_genomes'/f'{acc}.fna','fasta'))[contig]
        start=max(0,min(h['start'] for h in hits.values())-25000)
        end=min(len(seq),max(h['end'] for h in hits.values())+25000)
        rec=SeqRecord(seq.seq[start:end],id=f'{acc}_{contig}'[:50],name=acc[:16],
                      description=f'{acc} sol neighborhood {contig}:{start+1}-{end}')
        rec.annotations['molecule_type']='DNA'
        gene_by_protein={h['protein']:g for g,h in hits.items()}
        for pid,(c,s,e,strand,aa) in proteins(work,acc).items():
            if c!=contig or s<start or e>end:continue
            qualifiers={'locus_tag':[pid], 'product':['predicted protein'], 'translation':[aa]}
            if pid in gene_by_protein:
                qualifiers['gene']=[gene_by_protein[pid]]
                qualifiers['product']=['sol candidate homolog']
            rec.features.append(SeqFeature(FeatureLocation(s-start,e-start,strand=strand),type='CDS',qualifiers=qualifiers))
        out=work/'06_neighborhoods'/f'{acc}_sol_25kb.gb'
        SeqIO.write(rec,out,'genbank')
    print('Clinker-ready GenBank neighborhoods written to 06_neighborhoods/')


def flank_anchors(work):
    hits=locus_mapping(work)[PRIMARY]
    contig=next(iter(hits.values()))['contig']
    start=min(h['start'] for h in hits.values());end=max(h['end'] for h in hits.values())
    pmap=proteins(work,PRIMARY)
    left=sorted(((pid,v) for pid,v in pmap.items() if v[0]==contig and v[2]<=start and start-v[2]<=25000),
                key=lambda x:start-x[1][2])[:5]
    right=sorted(((pid,v) for pid,v in pmap.items() if v[0]==contig and v[1]>=end and v[1]-end<=25000),
                 key=lambda x:x[1][1]-end)[:5]
    if len(left)<2 or len(right)<2:raise RuntimeError('Fewer than 2 flank anchors on either side of primary sol locus')
    anchors={f'L{i}':p for i,p in enumerate(left,1)}
    anchors.update({f'R{i}':p for i,p in enumerate(right,1)})
    with (work/'09_flank_anchors'/'flank_queries.faa').open('w') as f:
        for name,(pid,v) in anchors.items():f.write(f'>{name} {pid}\n{v[4]}\n')
    write_tsv(work/'09_flank_anchors'/'primary_flank_anchors.tsv',
              [dict(anchor=n,protein=pid,contig=v[0],start=v[1],end=v[2]) for n,(pid,v) in anchors.items()],
              ['anchor','protein','contig','start','end'])
    rows=[]
    loci={r['accession']:r for r in read_tsv(work/'05_results'/'sol_loci.tsv')}
    for acc,row in loci.items():
        faa=work/'03_proteomes'/f'{acc}.faa';db=work/'04_blast'/f'{acc}_db'
        dest=work/'09_flank_anchors'/f'{acc}.tsv'
        if not dest.exists():
            with dest.open('w') as out:
                run(['blastp','-query',work/'09_flank_anchors'/'flank_queries.faa','-db',db,
                     '-evalue','1e-10','-max_target_seqs','10','-outfmt',
                     '6 qseqid sseqid pident qcovs evalue bitscore'],stdout=out)
        pmap=proteins(work,acc);best={}
        for line in dest.read_text().splitlines():
            a=line.split()
            if len(a)!=6:continue
            anchor,pid=a[:2];ident=float(a[2]);cov=float(a[3]);ev=float(a[4]);bits=float(a[5])
            if ident<35 or cov<60 or ev>1e-10 or pid not in pmap:continue
            if anchor not in best or bits>best[anchor][0]:best[anchor]=(bits,pid,pmap[pid])
        left_found=[(k,best[k][2]) for k in best if k.startswith('L')]
        right_found=[(k,best[k][2]) for k in best if k.startswith('R')]
        anchor_contigs=[v[2][0] for v in best.values()]
        same_contig=len(set(anchor_contigs))==1 if anchor_contigs else False
        gap='NA';orientation='NA'; verdict='not_assessable'
        if len(left_found)>=2 and len(right_found)>=2 and same_contig:
            # Check consistent order of L1/L2 and R1/R2; reversal allowed.
            L1=best.get('L1');R1=best.get('R1')
            if L1 and R1:
                ls,le=L1[2][1:3];rs,re=R1[2][1:3]
                gap=max(0,rs-le) if ls<rs else max(0,ls-re)
                orientation='forward' if ls<rs else 'reverse'
                verdict='candidate_empty_site' if row['status']=='no_coherent_cluster_detected' and gap<20000 else 'flanks_linked_inspect'
        rows.append(dict(accession=acc,sol_status=row['status'],left_anchors=len(left_found),right_anchors=len(right_found),
                         all_anchors_same_contig=same_contig,orientation=orientation,nearest_flank_gap_bp=gap,
                         preliminary_interpretation=verdict))
    write_tsv(work/'05_results'/'flank_anchor_comparison.tsv',rows,
              ['accession','sol_status','left_anchors','right_anchors','all_anchors_same_contig',
               'orientation','nearest_flank_gap_bp','preliminary_interpretation'])


def junctions(work):
    """Map two 2-kb sol-adjacent DNA anchors; report candidate empty sites, never HGT proof."""
    hits=locus_mapping(work)[PRIMARY]
    cid=next(iter(hits.values()))['contig']
    seq=str(SeqIO.to_dict(SeqIO.parse(work/'01_genomes'/f'{PRIMARY}.fna','fasta'))[cid].seq).upper()
    lo=min(h['start'] for h in hits.values());hi=max(h['end'] for h in hits.values())
    if lo<2000 or len(seq)-hi<2000:
        raise RuntimeError('Primary sol locus lacks 2-kb flanking sequence for DNA junction mapping')
    query=work/'09_flank_anchors'/'flank_DNA_2kb.fna'
    query.write_text(f'>LEFT_2kb\n{seq[lo-2000:lo]}\n>RIGHT_2kb\n{seq[hi:hi+2000]}\n')
    flanks={r['accession']:r for r in read_tsv(work/'05_results'/'flank_anchor_comparison.tsv')}
    loci={r['accession']:r for r in read_tsv(work/'05_results'/'sol_loci.tsv')}
    rows=[]
    for acc,row in loci.items():
        db=work/'09_flank_anchors'/f'{acc}_nucdb'
        output=work/'09_flank_anchors'/f'{acc}_DNA.tsv'
        if not output.exists():
            run(['makeblastdb','-in',work/'01_genomes'/f'{acc}.fna','-dbtype','nucl','-out',db],
                stdout=subprocess.DEVNULL)
            with output.open('w') as o:
                run(['blastn','-query',query,'-db',db,'-task','blastn','-evalue','1e-20',
                     '-max_target_seqs','10','-outfmt',
                     '6 qseqid sseqid pident qcovhsp sstart send bitscore'],stdout=o)
        best={}
        for line in output.read_text().splitlines():
            a=line.split()
            if len(a)!=7:continue
            q,contig=a[:2];ident=float(a[2]);cov=float(a[3]);s=int(a[4]);e=int(a[5]);bits=float(a[6])
            if ident<80 or cov<70:continue
            if q not in best or bits>best[q][0]:best[q]=(bits,contig,s,e,ident,cov)
        left=best.get('LEFT_2kb');right=best.get('RIGHT_2kb')
        same=bool(left and right and left[1]==right[1])
        gap='NA';orientation='NA';verdict='unresolved'
        if same:
            ls,le=left[2:4];rs,re=right[2:4]
            # Both anchors must map with the same orientation.
            if (ls<=le)==(rs<=re):
                if ls<rs:
                    gap=max(0,min(rs,re)-max(ls,le)-1);orientation='forward'
                else:
                    gap=max(0,min(ls,le)-max(rs,re)-1);orientation='reverse'
                f=flanks.get(acc,{})
                protein_support=(int(f.get('left_anchors',0))>=2 and int(f.get('right_anchors',0))>=2)
                if (row['status']=='no_coherent_cluster_detected' and gap<20000
                    and protein_support):verdict='candidate_empty_site_requires_alignment'
                elif row['status']=='candidate_intact_cluster':verdict='sol_positive_flank_reference'
                else:verdict='flanks_colocalized_inspect'
        rows.append(dict(accession=acc,sol_status=row['status'],left_DNA_match=bool(left),
                         right_DNA_match=bool(right),same_contig=same,orientation=orientation,
                         inter_anchor_gap_bp=gap,interpretation=verdict))
    write_tsv(work/'05_results'/'integration_site_candidates.tsv',rows,
              ['accession','sol_status','left_DNA_match','right_DNA_match','same_contig',
               'orientation','inter_anchor_gap_bp','interpretation'])


def composition(work):
    ref=json.loads((work/'02_reference'/'sol_reference.json').read_text())
    records=SeqIO.to_dict(SeqIO.parse(work/'01_genomes'/f'{PRIMARY}.fna','fasta'))
    cid=ref['contig'];
    if cid not in records:
        matches=[x for x in records if canonical_contig(x)==canonical_contig(cid)]
        if len(matches)!=1:raise RuntimeError('Reference contig not in primary genome')
        cid=matches[0]
    seq=str(records[cid].seq).upper()
    sol_nt=str(next(SeqIO.parse(work/'02_reference'/'sol_reference.fna','fasta')).seq).upper()
    from Bio.Seq import Seq
    reverse_nt=str(Seq(sol_nt).reverse_complement())
    forward_hits=seq.count(sol_nt);reverse_hits=seq.count(reverse_nt)
    if forward_hits+reverse_hits!=1:
        raise RuntimeError(f'Expected one exact sol DNA match in primary RefSeq, found {forward_hits} forward and {reverse_hits} reverse')
    start=seq.find(sol_nt) if forward_hits else seq.find(reverse_nt)
    end=start+len(sol_nt);size=end-start
    def gc(s):
        good=s.count('G')+s.count('C');at=s.count('A')+s.count('T')
        return 100*good/(good+at) if good+at else float('nan')
    left=seq[max(0,start-25000):start];right=seq[end:end+25000]
    rng=random.Random(20261008)
    pool=[(rid,str(r.seq).upper()) for rid,r in records.items() if len(r)>=size]
    weights=[len(s)-size+1 for rid,s in pool]
    draws=[];attempts=0
    while len(draws)<10000 and attempts<200000:
        attempts+=1;idx=rng.choices(range(len(pool)),weights=weights,k=1)[0]
        rid,s=pool[idx];pos=rng.randrange(len(s)-size+1)
        if rid==cid and pos<end and pos+size>start:continue
        draws.append(gc(s[pos:pos+size]))
    if len(draws)<1000:raise RuntimeError('Insufficient genomic background windows for GC comparison')
    obs=gc(seq[start:end]);avg=sum(draws)/len(draws)
    sd=math.sqrt(sum((x-avg)**2 for x in draws)/(len(draws)-1))
    p=(1+sum(abs(x-avg)>=abs(obs-avg) for x in draws))/(len(draws)+1)
    data={'sol_gc_pct':round(obs,3),'genome_gc_pct':round(gc(''.join(str(r.seq).upper() for r in records.values())),3),
          'left_25kb_gc_pct':round(gc(left),3),'right_25kb_gc_pct':round(gc(right),3),
          'random_window_mean_gc_pct':round(avg,3),'random_window_sd_gc_pct':round(sd,3),
          'sol_gc_z_vs_windows':round((obs-avg)/sd,3) if sd else 'NA',
          'two_sided_empirical_p':round(p,6),'n_background_windows':len(draws),
          'note':'Exploratory, non-independent genomic windows; composition alone does not prove HGT.'}
    (work/'08_composition'/'GC_permutation.json').write_text(json.dumps(data,indent=2)+'\n')
    print('GC composition:',data)


def mobile_features(work):
    # Annotated GBFF, not Prodigal labels; unavailable genomes are explicitly NA.
    pattern=re.compile(r'integrase|transposase|insertion sequence|recombinase|phage|holin|lysozyme|tRNA',re.I)
    loci={r['accession']:r for r in read_tsv(work/'05_results'/'sol_loci.tsv')}
    rows=[]
    for acc,row in loci.items():
        if row['status']!='candidate_intact_cluster':continue
        gbff=work/'01_genomes'/f'{acc}.gbff'
        if acc==PRIMARY:
            candidates=list((work/'02_reference'/'source'/'package'/'ncbi_dataset'/'data').rglob('*.gbff'))
            if candidates:gbff=candidates[0]
        if not gbff.exists():
            rows.append(dict(accession=acc,annotation_available='no',feature='NA',product='NA',distance_to_sol_bp='NA'))
            continue
        target=row['contig'];lo=int(row['start_1based'])-1;hi=int(row['end']);count=0
        for rec in SeqIO.parse(gbff,'genbank'):
            if canonical_contig(rec.id)!=canonical_contig(target):continue
            for feat in rec.features:
                if feat.type not in ('CDS','tRNA','tmRNA'):continue
                s=int(feat.location.start);e=int(feat.location.end)
                if e<lo-25000 or s>hi+25000:continue
                prod='; '.join(feat.qualifiers.get('product',[]))
                if feat.type=='tRNA' or pattern.search(prod):
                    distance=max(lo-e,s-hi,0)
                    rows.append(dict(accession=acc,annotation_available='yes',feature=feat.type,
                                     product=prod or 'tRNA',distance_to_sol_bp=distance))
                    count+=1
        if count==0:rows.append(dict(accession=acc,annotation_available='yes',feature='none_detected',
                                     product='NA',distance_to_sol_bp='NA'))
    write_tsv(work/'05_results'/'mobile_element_candidates.tsv',rows,
              ['accession','annotation_available','feature','product','distance_to_sol_bp'])


def make_report(work):
    rows=read_tsv(work/'05_results'/'sol_loci.tsv')
    positives=[r for r in rows if r['status']=='candidate_intact_cluster']
    genus=[r for r in rows if r['organism'].startswith('Rouxiella ')]
    text=(f'# Preliminary sol HGT analysis (not a conclusion)\n\n'
          f'- Genomes examined: {len(rows)}; Rouxiella: {len(genus)}\n'
          f'- Candidate compact sol loci: {len(positives)}\n'
          '- Absence in fragmented assemblies is provisional.\n'
          '- Homology calls require confirmation by reciprocal similarity, gene order and domain architecture.\n'
          '- A GC shift is supporting evidence only; verify integration junctions before asserting HGT.\n'
          '- Published Matilla et al. (2022) already reported intergeneric sol HGT evidence.\n'
          '- Compare trees only on identical taxa; discordance alone is not proof of HGT.\n\n'
          '## Files\n'
          '- `05_results/sol_loci.tsv`, `sol_gene_matrix.tsv`, `sol_gene_hits.tsv`\n'
          '- `05_results/flank_anchor_comparison.tsv`, `integration_site_candidates.tsv`, `mobile_element_candidates.tsv`\n'
          '- `06_neighborhoods/*.gb` (Clinker-ready)\n'
          '- `08_composition/GC_permutation.json`\n')
    (work/'05_results'/'README_results.md').write_text(text)


def discover(root,work):
    complete(work,'01_genomes',lambda:discover_genomes(root,work),[work/'05_results'/'genome_manifest.tsv'])
    complete(work,'02_reference',lambda:extract_reference(work),[work/'02_reference'/'sol_queries.faa',work/'02_reference'/'sol_reference.json'])
    complete(work,'03_proteomes',lambda:predict(work),[work/'03_proteomes'/f'{PRIMARY}.faa'])
    complete(work,'04_blast',lambda:blast(work),[work/'04_blast'/f'{PRIMARY}.tsv'])
    complete(work,'05_summary',lambda:summarize(work),[work/'05_results'/'sol_loci.tsv',work/'05_results'/'sol_gene_hits.tsv'])
    complete(work,'06_neighborhoods',lambda:neighborhoods(work),[work/'06_neighborhoods'/f'{PRIMARY}_sol_25kb.gb'])
    complete(work,'07_flank_anchors',lambda:flank_anchors(work),[work/'05_results'/'flank_anchor_comparison.tsv'])
    complete(work,'08_junctions',lambda:junctions(work),[work/'05_results'/'integration_site_candidates.tsv'])
    complete(work,'09_composition',lambda:composition(work),[work/'08_composition'/'GC_permutation.json'])
    complete(work,'10_mobile_elements',lambda:mobile_features(work),[work/'05_results'/'mobile_element_candidates.tsv'])
    complete(work,'11_report',lambda:make_report(work),[work/'05_results'/'README_results.md'])


def align_sequences(src,dst):
    with dst.open('w') as o:run(['mafft','--auto',src],stdout=o)


def trim(src,dst):
    with dst.open('w') as o:run(['trimal','-in',src,'-automated1'],stdout=o)
    if dst.stat().st_size==0:raise RuntimeError(f'Alignment empty after trimAl: {src}')


def phylo_orthofinder(work,threads):
    pdir=work/'07_phylogeny'/'proteomes';pdir.mkdir(parents=True,exist_ok=True)
    for r in read_tsv(work/'05_results'/'genome_manifest.tsv'):
        acc=r['accession']; dest=pdir/f'{acc}.fa';src=work/'03_proteomes'/f'{acc}.faa'
        if not dest.exists():dest.symlink_to(src.resolve())
    run(['orthofinder','-f',pdir,'-t',str(threads),'-a',str(threads),'-M','msa'])
    result=sorted((pdir/'OrthoFinder').glob('Results_*'))
    if not result:raise RuntimeError('OrthoFinder results not found')
    (work/'07_phylogeny'/'orthofinder_results_path.txt').write_text(str(result[-1])+'\n')


def phylo_species(work,threads):
    p=Path((work/'07_phylogeny'/'orthofinder_results_path.txt').read_text().strip())
    ortho=p/'Orthogroups'/'Orthogroups.tsv'
    if not ortho.exists():raise RuntimeError(f'No OrthoFinder orthogroup table: {ortho}')
    with ortho.open() as f:
        reader=csv.DictReader(f,delimiter='\t');species=[x for x in reader.fieldnames if x!='Orthogroup']
        single=[]
        for r in reader:
            if all(r[s].strip() and ',' not in r[s] for s in species):single.append(r)
    if len(single)<100:raise RuntimeError(f'Only {len(single)} single-copy orthogroups; inspect OrthoFinder')
    # Deterministic selection of up to 200 single-copy orthogroups.
    random.Random(20261008).shuffle(single);single=single[:min(200,len(single))]
    sequences={s:SeqIO.to_dict(SeqIO.parse(work/'03_proteomes'/f'{s}.faa','fasta')) for s in species}
    aln_dir=work/'07_phylogeny'/'core_alignments';aln_dir.mkdir(exist_ok=True)
    alignments=[]
    for r in single:
        og=r['Orthogroup'];f=aln_dir/f'{og}.fa'
        with f.open('w') as out:
            for s in species:
                pid=r[s].strip();rec=sequences[s].get(pid)
                if rec is None:raise RuntimeError(f'OrthoFinder protein {pid} absent from {s}')
                out.write(f'>{s}\n{rec.seq}\n')
        aligned=aln_dir/f'{og}.aligned.fa'; trimmed=aln_dir/f'{og}.trimmed.fa'
        align_sequences(f,aligned);trim(aligned,trimmed)
        alignments.append(trimmed)
    concat={s:'' for s in species};partitions=[];pos=1
    for f in alignments:
        records=SeqIO.to_dict(SeqIO.parse(f,'fasta'))
        length=len(next(iter(records.values())).seq)
        if length<10:continue
        if set(records)!=set(species):raise RuntimeError(f'Missing species in alignment {f}')
        for s in species:concat[s]+=str(records[s].seq)
        partitions.append((f.stem.split('.')[0],pos,pos+length-1));pos+=length
    target=work/'07_phylogeny'/'species_core_concat.fa'
    with target.open('w') as out:
        for s,seq in concat.items():out.write(f'>{s}\n{seq}\n')
    write_tsv(work/'07_phylogeny'/'species_partitions.tsv',
              [dict(gene=x,start=y,end=z) for x,y,z in partitions],['gene','start','end'])
    print(f'Core alignment: {len(partitions)} orthogroups, {pos-1} aa sites')
    if pos-1<10000:raise RuntimeError('Too few conserved positions for species tree')
    run(['iqtree2','-s',target,'-m','LG+F+G4','-B','1000','-T',str(threads),
         '--prefix',work/'07_phylogeny'/'species_core'])


def phylo_sol(work,threads):
    mapping=locus_mapping(work)
    rows={r['accession']:r for r in read_tsv(work/'05_results'/'sol_loci.tsv')}
    eligible=[a for a,h in mapping.items() if rows[a]['status']=='candidate_intact_cluster' and
              sum(g in h for g in SOL_TREE_GENES)>=5]
    genus={rows[a]['organism'].split()[0] for a in eligible}
    if len(eligible)<5 or len(genus)<2:
        raise RuntimeError(f'Insufficient sol-positive genomes/taxonomic diversity for tree: {len(eligible)}, genera={genus}')
    seqmaps={a:proteins(work,a) for a in eligible}
    outdir=work/'07_phylogeny'/'sol_alignments';outdir.mkdir(exist_ok=True)
    concat={a:'' for a in eligible};parts=[];pos=1
    for gene in SOL_TREE_GENES:
        f=outdir/f'{gene}.fa';present=[a for a in eligible if gene in mapping[a]]
        if len(present)<4:continue
        with f.open('w') as out:
            for a in present:
                pid=mapping[a][gene]['protein'];out.write(f'>{a}\n{seqmaps[a][pid][4]}\n')
        aligned=outdir/f'{gene}.aligned.fa';trimmed=outdir/f'{gene}.trimmed.fa'
        align_sequences(f,aligned);trim(aligned,trimmed)
        seqs=SeqIO.to_dict(SeqIO.parse(trimmed,'fasta'))
        n=len(next(iter(seqs.values())).seq)
        if n<20:continue
        for a in eligible:concat[a]+=str(seqs[a].seq) if a in seqs else '-'*n
        parts.append((gene,pos,pos+n-1));pos+=n
    if len(parts)<4:raise RuntimeError('Too few usable orthologous sol genes; inspect alignments')
    target=work/'07_phylogeny'/'sol_concat.fa'
    with target.open('w') as out:
        for a in eligible:out.write(f'>{a}\n{concat[a]}\n')
    write_tsv(work/'07_phylogeny'/'sol_partitions.tsv',
              [dict(gene=g,start=s,end=e) for g,s,e in parts],['gene','start','end'])
    (work/'07_phylogeny'/'sol_tree_taxa.txt').write_text('\n'.join(eligible)+'\n')
    print(f'Sol alignment: {len(parts)} genes, {pos-1} aa sites, {len(eligible)} taxa')
    run(['iqtree2','-s',target,'-m','LG+F+G4','-B','1000','-T',str(threads),
         '--prefix',work/'07_phylogeny'/'sol_cluster'])


def phylo_gene_trees(work,threads):
    """Independent gene trees provide a concordance check on the concatenated sol topology."""
    d=work/'07_phylogeny'/'sol_alignments'
    rows=[]
    for gene in SOL_TREE_GENES:
        aligned=d/f'{gene}.trimmed.fa'
        if not aligned.exists():continue
        records=list(SeqIO.parse(aligned,'fasta'))
        length=len(records[0].seq) if records else 0
        if len(records)<5 or length<60:
            rows.append(dict(gene=gene,taxa=len(records),sites=length,status='too_little_information'))
            continue
        prefix=d/f'{gene}_ML'
        if not Path(str(prefix)+'.treefile').exists():
            run(['iqtree2','-s',aligned,'-m','LG+F+G4','-B','1000','-T',str(threads),'--prefix',prefix])
        rows.append(dict(gene=gene,taxa=len(records),sites=length,status='tree_inferred'))
    write_tsv(work/'07_phylogeny'/'individual_sol_gene_trees.tsv',rows,
              ['gene','taxa','sites','status'])


def splits(tree):
    names=set(t.name for t in tree.get_terminals());n=len(names);out=set()
    for node in tree.get_nonterminals():
        a=frozenset(t.name for t in node.get_terminals());b=frozenset(names-a)
        if len(a)<2 or len(b)<2:continue
        out.add(min((a,b),key=lambda z:(len(z),tuple(sorted(z)))))
    return out


def compare_trees(work,threads):
    d=work/'07_phylogeny'
    species=Phylo.read(d/'species_core.treefile','newick')
    sol=Phylo.read(d/'sol_cluster.treefile','newick')
    taxa=set(t.name for t in sol.get_terminals())
    for t in list(species.get_terminals()):
        if t.name not in taxa:species.prune(t)
    if set(t.name for t in species.get_terminals())!=taxa:raise RuntimeError('Species and sol tree taxa differ')
    Phylo.write(species,d/'species_pruned_to_sol_taxa.nwk','newick')
    a=splits(species);b=splits(sol);rf=len(a.symmetric_difference(b))
    result={'taxa':len(taxa),'species_splits':len(a),'sol_splits':len(b),
            'shared_splits':len(a&b),'unrooted_robinson_foulds_distance':rf,
            'interpretation':'RF discordance alone is not proof of HGT; check branch support and AU test.'}
    (d/'tree_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    with (d/'topology_candidates.nwk').open('w') as f:
        f.write((d/'sol_cluster.treefile').read_text().strip()+'\n')
        f.write((d/'species_pruned_to_sol_taxa.nwk').read_text().strip()+'\n')
    run(['iqtree2','-s',d/'sol_concat.fa','-m','LG+F+G4','-z',d/'topology_candidates.nwk',
         '-n','0','-zb','1000','-au','-T',str(threads),'--prefix',d/'topology_AU_test'])
    print('Tree comparison:',result)


def phylogeny(work,threads):
    if not (work/'05_results'/'sol_gene_hits.tsv').exists():
        raise RuntimeError('Run discover stage first')
    complete(work,'12_orthofinder',lambda:phylo_orthofinder(work,threads),
             [work/'07_phylogeny'/'orthofinder_results_path.txt'])
    complete(work,'13_species_tree',lambda:phylo_species(work,threads),
             [work/'07_phylogeny'/'species_core.treefile'])
    complete(work,'14_sol_tree',lambda:phylo_sol(work,threads),
             [work/'07_phylogeny'/'sol_cluster.treefile'])
    complete(work,'15_gene_trees',lambda:phylo_gene_trees(work,threads),
             [work/'07_phylogeny'/'individual_sol_gene_trees.tsv'])
    complete(work,'16_tree_comparison',lambda:compare_trees(work,threads),
             [work/'07_phylogeny'/'tree_comparison.json',work/'07_phylogeny'/'topology_AU_test.iqtree'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['discover','phylogeny'])
    p.add_argument('--root',type=Path,default=Path('/scratch/al98750/Roux'))
    p.add_argument('--threads',type=int,default=8)
    a=p.parse_args();w=setup(a.root)
    if a.stage=='discover':discover(a.root,w)
    else:phylogeny(w,a.threads)
    print(f'Stage {a.stage} complete; outputs: {w}')

if __name__=='__main__':main()
