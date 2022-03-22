#!/bin/bash

#Para correr este código se necesita: 
# - genoma del virus en fasta (covid-wuhan.fasta)
# - índices del genoma del virus (incluído en la carpeta "programas")
# - indices Nextera (si se hace el secuenciamiento std Illumina [NexteraPE-PE.fa])
# - Trimmomatic (Incluído en la carpeta "programas")

############# Definition of N° cores ###########

threads=8

############## Definición de muestras a trabajar ####### 

#Lista de muestras forward y reverse
ls secuencias/*R1_001.fastq.gz > list.txt 

#Me quedo con el nombre de la muestra, quitando los prefijos que indican si la muestra es reverse o forward
list=$( sed -e 's/_R1_001.fastq.gz//' -e 's/secuencias//' -e 's/[/]//g' list.txt )

#conversión de la lista a un array
arr=($list)

#Se cuentan el número de muestras con las que se va a trabajar
howmany() { echo $#; }

a=$( howmany $list )

b=$(expr $a - 1)

#Creación de carpetas donde van a ir todos los resultados de la corrida 

mkdir Resultados  

cp Programs/truseqht-PE.fa Resultados/

cd Resultados
#Aquí se almacenan los fastq que hicieron match con el genoma de referencia de Sars-cov-2 de wuhan

mkdir solo_match_sars-cov-2

#Aquí se almacenan todos los genomas ensamblados de megahit (fasta + muchos productos más generados por el programa). 
mkdir ensamble-megahit

#Aquí se almacenan todos los contigs ensamblados en fasta.
mkdir contigs-megahit

#Aquí se almacenan los reportes de coverage de las secuencias 
mkdir coverage 

############################ Protocolo ###########################

for ((c=0;c<=$b;c++))
do 

#Filtro de calidad de las muestras 
java -jar ../Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE\
   -threads $threads \
   ../secuencias/${arr[$c]}_R1_001.fastq.gz \
   ../secuencias/${arr[$c]}_R2_001.fastq.gz \
   ${arr[$c]}-R1_paired.fastq \
   ${arr[$c]}-R1_unpaired.fastq \
   ${arr[$c]}-R2_paired.fastq \
   ${arr[$c]}-R2_unpaired.fastq \
   ILLUMINACLIP:truseqht-PE.fa:2:30:10:2:keepBothReads LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30  

#Enfrentamiento entre la muestra con el genoma de referencia de Sars-cov-2 (Muestra de wuhan). Nos quedamos con aquellas secuencias que hayan hecho match con nuestras muestras 

../Programs/bowtie2-2.2.6/bowtie2 -x ../Programs/COVID-reference/virus \
  -1 ${arr[$c]}-R1_paired.fastq \
  -2 ${arr[$c]}-R2_paired.fastq \
  -U ${arr[$c]}-R1_unpaired.fastq,${arr[$c]}-R2_unpaired.fastq \
  -p $threads \
  --no-unal \
  -S  ${arr[$c]}-SAMPLE_mapped.sam 

#Conversión de sam a bam de las secuencias que hicieron match 
../Programs/samtools-1.9/samtools view \
  -bS ${arr[$c]}-SAMPLE_mapped.sam \
  --threads $threads \
  > ${arr[$c]}-SAMPLE_mapped.bam

#Orden de las muestras
../Programs/samtools-1.9/samtools sort \
  ${arr[$c]}-SAMPLE_mapped.bam \
  --threads $threads \
  -o ${arr[$c]}-SAMPLE_mapped_sorted.bam

../Programs/samtools-1.9/samtools fastq \
-0 /dev/null \
-s ${arr[$c]}-single.fq \
-N ${arr[$c]}-SAMPLE_mapped_sorted.bam \
--threads $threads \
> ${arr[$c]}-paired.fq 

mv *.fq solo_match_sars-cov-2 

cd ensamble-megahit

#Ensamblaje de las secuencias filtradas 

../../Programs/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit \
--12 ../solo_match_sars-cov-2/${arr[$c]}-paired.fq \
-r ../solo_match_sars-cov-2/${arr[$c]}-single.fq \
-o ${arr[$c]}-ensamble \
-t $threads 

cd ..
done

#Copiar los contigs de los ensamblajes de megahit a la carpeta "contigs-megahit" 

for ((c=0;c<=$b;c++))
do 
cp -r ensamble-megahit/${arr[$c]}-ensamble/final.contigs.fa contigs-megahit/
cd contigs-megahit 
  mv final.contigs.fa  ${arr[$c]}.fasta 
cd ..
done 

#Conteo de genomas que están (o no) completos (y que se tienen que salvar) 

echo '' > genomas_completos.txt

for ((c=0;c<=$b;c++))
do 

d=$( grep '>' contigs-megahit/${arr[$c]}.fasta | wc -l )

if [ $d != 1 ] 
then
	echo 'El genoma '${arr[$c]}.fasta' no está completo' >> genomas_completos.txt
else 
	echo 'El genoma '${arr[$c]}.fasta' está completo' >> genomas_completos.txt
fi 

done

#Cálculo del coverage 

mv *_sorted.bam coverage/ 

cd coverage 

for ((c=0;c<=$b;c++))
do 

if [ $c == 0 ] 
then 
	../../Programs/samtools-1.9/samtools depth -a ${arr[$c]}-SAMPLE_mapped_sorted.bam |  awk '{sum+=$3} END { print "'${arr[$c]}' coverage = ",sum/29903}' > coverage.txt 
else 
	../../Programs/samtools-1.9/samtools depth -a ${arr[$c]}-SAMPLE_mapped_sorted.bam |  awk '{sum+=$3} END { print "'${arr[$c]}' coverage = ",sum/29903}' >> coverage.txt 
fi
done

