out=${1#"/fastdata/bop20pp/supergene_AB_maleref/gvcfs/NC_044241.2"}
out=${out%".g.vcf.gz"}
echo $out

echo $out > /fastdata/bop20pp/supergene_AB_maleref/gvcfs/combine/${out}_reheader.txt

bcftools reheader \
	-s /fastdata/bop20pp/supergene_AB_maleref/gvcfs/combine/${out}_reheader.txt \
	-o /fastdata/bop20pp/supergene_AB_maleref/gvcfs/NC_044241.2${out}.rh.g.vcf.gz $1
