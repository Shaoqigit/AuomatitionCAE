Introduction
======
This is a project accomplished in the Institut de Mécanique et d'Ingéniere de Bordeaux (I2M) concerning the theorical analysis and numerical calculation of adhesively bonded joints.  
Analytical implementation uses the volkersen model and numercail analysis is studied in Abaqus with different methods.  
All of these code is targeted to evaluate the faillure probability of any mechanical structure basing methode AK-MCS.  
# Softwares and language used
Python script/spyder  
CATIA v5 6R2018  
Abaqus v2016  
# Instruction
1. The file named stress_analyse is based on the Volkersen model to analyse stress in the cohsive zone. Two python script, one is about using analytical performance fucntion of Volkersen model in evaluation process, while the performace function is computered by abaqus in another script.  
![image](https://github.com/Shaoqigit/Automatisation-tool-for-evaluating-of-adhesively-bonded/blob/master/figure/1.png)
2. The file named CATIA_with_Python can generate NACA profil series 4 in CATIA with equations coded in Python
![image](https://github.com/Shaoqigit/Automatisation-tool-for-evaluating-of-adhesively-bonded/blob/master/CATIA_with_Python/profil_catia.PNG)  
Another file is to extract all coordiantes in catia saved in Excel
3. The file named Cohsive_zone_model can analyse adhesivedly bonded joints with <Cohesive Zone Models> in Abaqus.  
  ![image](https://github.com/Shaoqigit/Automatisation-tool-for-evaluating-of-adhesively-bonded/blob/master/Cohsive_zone_model/CZMs.png) 
