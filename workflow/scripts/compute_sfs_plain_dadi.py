import dadi
import fastdfe as fd
import tempfile

try:
    testing = False
    vcf = snakemake.input.vcf
    ingroups = snakemake.params["ingroups"]
    n = snakemake.params["n"]
    out = snakemake.output[0]
except NameError:
    testing = True
    vcf = "results/vcf/Homo_sapiens/Pongo_abelii.degeneracy.polarized.vcf.gz"
    ingroups = ["gorilla_gorilla_gorilla"]
    n = 20
    out = "scratch/chr22.sfs.20.dadi.fs"

# Write temporary popinfo file
with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
    for name in ingroups:
        tmp.write(f"{name}\tpop0\n")
    tmp.flush()

    dd = dadi.Misc.make_data_dict_vcf(vcf, popinfo_filename=tmp.name)

fs = dadi.Spectrum.from_data_dict(dd, pop_ids=["pop0"], projections=[n], polarized=True)

sfs = fd.Spectrum(fs.data)
sfs.to_file(out)