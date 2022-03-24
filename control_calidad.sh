#!/bin/bash

#Para correr este código se necesita: 
# - genoma del virus en fasta (covid-wuhan.fasta)
# - índices del genoma del virus (incluído en la carpeta "Programas") 
# - indices Truseq (Esto depende de los índices utilizados en tu secuenciamiento. En caso uses otros adaptadores, colocar el formato fasta en la carpeta programas y cambiar "truseqht-PE" por el nombre de tu adaptador en la línea 65 del código. 

############# Variables que puedes modificar ###########

threads=4

adapters="truseqht-PE.fa"

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

#Selecciona el fasta de tus adaptadores y mover a la carpeta de trabajo
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
trimmomatic PE\
   -threads $threads \
   ../secuencias/${arr[$c]}_R1_001.fastq.gz \
   ../secuencias/${arr[$c]}_R2_001.fastq.gz \
   ${arr[$c]}-R1_paired.fastq \
   ${arr[$c]}-R1_unpaired.fastq \
   ${arr[$c]}-R2_paired.fastq \
   ${arr[$c]}-R2_unpaired.fastq \
   ILLUMINACLIP:$adapters:2:30:10:2:keepBothReads LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30  

#Enfrentamiento entre la muestra con el genoma de referencia de Sars-cov-2 (Muestra de wuhan). Nos quedamos con aquellas secuencias que hayan hecho match con nuestras muestras 

bowtie2 -x ../Programs/COVID-reference/virus \
  -1 ${arr[$c]}-R1_paired.fastq \
  -2 ${arr[$c]}-R2_paired.fastq \
  -U ${arr[$c]}-R1_unpaired.fastq,${arr[$c]}-R2_unpaired.fastq \
  -p $threads \
  --no-unal \
  -S  ${arr[$c]}-SAMPLE_mapped.sam 

#Conversión de sam a bam de las secuencias que hicieron match 
samtools view \
  -bS ${arr[$c]}-SAMPLE_mapped.sam \
  --threads $threads \
  > ${arr[$c]}-SAMPLE_mapped.bam

#Orden de las muestras
samtools sort \
  ${arr[$c]}-SAMPLE_mapped.bam \
  --threads $threads \
  -o ${arr[$c]}-SAMPLE_mapped_sorted.bam

samtools fastq \
-0 /dev/null \
-s ${arr[$c]}-single.fq \
-N ${arr[$c]}-SAMPLE_mapped_sorted.bam \
--threads $threads \
> ${arr[$c]}-paired.fq 

mv *.fq solo_match_sars-cov-2 

cd ensamble-megahit

#Ensamblaje de las secuencias filtradas 

megahit \
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

cd contigs-megahit/

mkdir genomas_completos genomas_incompletos

cd ../

#Asignación de genomas completos

for ((c=0;c<=$b;c++))
do 

d=$( grep '>' contigs-megahit/${arr[$c]}.fasta | wc -l )

if [ $d != 1 ] 
then
	echo 'El genoma '${arr[$c]}.fasta' no está completo' >> genomas_completos.txt
	
	cd contigs-megahit/

	mv ${arr[$c]}.fasta genomas_incompletos/

	cd ../

else 
	echo 'El genoma '${arr[$c]}.fasta' está completo' >> genomas_completos.txt

	cd contigs-megahit/

	mv ${arr[$c]}.fasta genomas_completos/

	cd ../
fi 

done

#Cálculo del coverage 

mv *_sorted.bam coverage/ 

cd coverage 

for ((c=0;c<=$b;c++))
do 

if [ $c == 0 ] 
then 
	samtools depth -a ${arr[$c]}-SAMPLE_mapped_sorted.bam |  awk '{sum+=$3} END { print "'${arr[$c]}' coverage = ",sum/29903}' > coverage.txt 
else 
	samtools depth -a ${arr[$c]}-SAMPLE_mapped_sorted.bam |  awk '{sum+=$3} END { print "'${arr[$c]}' coverage = ",sum/29903}' >> coverage.txt 
fi
done

