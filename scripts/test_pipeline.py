import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from io import StringIO

p=Path(__file__).parent/'03_sol_HGT.py'
spec=importlib.util.spec_from_file_location('solhgt',p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
assert len(m.SOL_TAGS)==13
assert m.SOL_TAGS['solM']=='LLR01_11545'
assert m.canonical_contig('NZ_JAJGAU010000005.1')==m.canonical_contig('JAJGAU010000005.1')
# Prevent double counting one generic multidomain protein as 2 sol genes.
hits=[dict(gene='solA',protein='p1',contig='ctg',start=100,end=200,bitscore=100),
      dict(gene='solB',protein='p1',contig='ctg',start=100,end=200,bitscore=99),
      dict(gene='solB',protein='p2',contig='ctg',start=220,end=300,bitscore=90)]
contig, chosen=m.select_locus(hits)
assert contig=='ctg' and set(chosen)=={'solA','solB'} and chosen['solB']['protein']=='p2'
t1=Phylo.read(StringIO('((A,B),(C,(D,E)));'),'newick')
t2=Phylo.read(StringIO('((A,C),(B,(D,E)));'),'newick')
assert len(m.splits(t1).symmetric_difference(m.splits(t2)))>0
with TemporaryDirectory() as tmp:
    root=Path(tmp);work=m.setup(root)
    import random
    rng=random.Random(42)
    seq=''.join(rng.choices('ACGT',k=90000))
    rec=SeqRecord(Seq(seq),id='JAJGAU010000005.1',description='test')
    rec.annotations['molecule_type']='DNA'
    for i,(gene,tag) in enumerate(m.SOL_TAGS.items()):
        start=2000+i*900;end=start+300
        rec.features.append(SeqFeature(FeatureLocation(start,end,strand=1),type='CDS',
                            qualifiers={'locus_tag':[tag],'translation':['M'*99]}))
    pack=work/'tmp'/'source_package';pack.mkdir(parents=True)
    gbff=pack/'source.gbff';SeqIO.write(rec,gbff,'genbank')
    old=m.download
    m.download=lambda *args,**kwargs:pack
    m.extract_reference(work)
    m.download=old
    ref=m.json.loads((work/'02_reference'/'sol_reference.json').read_text())
    assert ref['contig']=='JAJGAU010000005.1' and len(ref['genes'])==13
    genome=work/'01_genomes'/f'{m.PRIMARY}.fna'
    SeqIO.write(SeqRecord(Seq(seq),id='NZ_JAJGAU010000005.1',description='test'),genome,'fasta')
    m.composition(work)
    gc=m.json.loads((work/'08_composition'/'GC_permutation.json').read_text())
    assert gc['n_background_windows']==10000
print('PASS: 13 sol tags, contig normalization, nonduplicated locus assignment, tree splits, GenBank reference extraction, GC window analysis')
