Ensamblaje de secuencias SARS-CoV-2

- Este script está diseñados para ensamblar secuencias crudas paired-end de SARS-CovV-2 secuenciadas con tecnología Illumina.

- Necesitas crear un ambiente en conda, donde se instalarán todos los programas necesarios para ensamblar genomas de SARS-CoV-2. 

- El script asume que tienes un procesador de 4 núcleos. Puedes modificar este valor entrando al archivo "control_calidad.sh" y cambiando el valor de la variable "threads". También asume que usaste adaptadores truseq para tu secuenciamiento. Puedes cambiar estos adaptadores por los de tu interés agregándolos en la carpeta "Programs" y escribiendo su nombre exacto en la variable "adapters"

Instalación

1. Instalación de conda (en caso no lo tengas instalado)

> wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

> chmod +x Miniconda3-latest-Linux-x86_64.sh

> ./Miniconda3-latest-Linux-x86_64.sh

2. Instalación del programa

> git clone https://github.com/leonardo467/Enasamblaje-SARS-CoV-2.git

> conda env create -n covid -f environment_covid.yaml

> rm environment_covid.yaml

3. Haz ejecutable el script

> chmod +x control_calidad.sh

Instrucciones

1. Copiar tus secuencias de interés (paired-end) en una carpeta llamada "secuencias". Es importante que el nombre sea "secuencias" para que pueda correr el código. Puedes descargar esta carpeta de secuencias de prueba para probar el script:

https://drive.google.com/drive/folders/1-8ZCpA83cyRR6ITGU8HVLy1itwheFkPr?usp=sharing

2. Activa el ambiente conda:

> conda activate covid

4. Ejecuta el script: 

> ./control_calidad.sh

*Mayor información y tutorial de cómo correr este programa se encuentra en el manual. 

Cualquier consulta/sugerencia sobre este script, pueden comunicarse a los siguientes correos:

Luis González --> luis.gonzalez.v@upch.pe

Diego Cuicapuza --> diego.cuicapuza@upch.pe

