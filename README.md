# Enasamblaje-SARS-CoV-2

############################# notas ############################################

- Este script está diseñados para ensamblar secuencias crudas paired-end de SARS-CovV-2 paired-end (Secuenciadas con tecnología Illumina) 

- Necesitas instalar (o verificar que tienes instalado) unos programas antes de comenzar. La lista de programas está al final de este texto. 

- El script asume que tienes un procesador de 4 núcleos. Puedes modificar este valor entrando al archivo "control_calidad.sh" y cambiando el valor de la variable 

################# Instrucciones de uso del script ###############################

 1. Copiar tus secuencias de interés (paired-end) en una carpeta llamada "secuencias". Es importante que el nombre sea "secuencias" para que pueda correr el código. 

 2. Activar la carpeta de programas para que todos los archivos sean ejecutables (esto solo se hace la primera vez que se usa el script). También es necesario hacer ejecutable el programa "control_calidad.sh" (Esto también se hace solo una vez en la vida). Este se hace de la siguiente manera: 

	> chmod -R +x Programs/
	
	> chmod +x control_calidad.sh
	
 3. Ejecutar el programa de la siguiente manera: 
 
 	> ./control_calidad.sh 

################# Programas que necesitas instalar ##############################

La gran mayoría de ellos ya vienen en la carpeta "Programs". Sin embargo, hay algunos que necesitas instalar para que pueda funcionar el script: 

- Java        --> Código de instalación --> sudo apt install openjdk-14-jre-headless 

- libncurses5 --> Código de intalación --> sudo apt-get install libncurses5

*Cualquier consulta/sugerencia sobre este script, pueden comunicarse al siguiente correo: 

       > Luis González --> luis.gonzalez.v@upch.pe
