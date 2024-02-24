## Global Genomic Network Communities (GG-NC) Browser

This GitHub page contains the code, input data and the wrapper to run the Shiny app for the [Global Genomic Network Communities Browser](https://sohail-lab.shinyapps.io/GG-NC/) developed by the [Sohail Lab](https://www.sohaillab.com/) at the Center for Genomic Sciences - UNAM, financed by Conahcyt.

Input data are obtained from the [1000 Genomes Project](https://www.nature.com/articles/nature15393) and [Human Genome Diversity Project](https://www.science.org/doi/10.1126/science.aay5012?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed).

The Shiny app, was developed within the research titled "Beyond Continental Groups for Framing Human Diversity: shifting from Static Clusters to Dynamic Communities in Genetic Networks", aims to make the results of this research more interactive and allow engagement from the community of scientists and the general public alike. This Shiny app allows you to view the different communities that emerge at different resolution values and their geographic distribution considering different genetic metrics (Genetic Relationship Matrix based on common or rare variants, Identity By Descent and Principal Component Analysis)  

## Shiny interface
Follow [this link](https://sohail-lab.shinyapps.io/GG-NC/) for the interactive Shiny app. A screenshot of the interface is provided below.

![Shiny app interface](app_image.jpg)

The Shiny app contains a toolbar with the following options:

+ __Select a genetic metric__: You can choose one of the provided genetic metrics: Genetic Relationship Matrix (GRM) for rare or common variants, Identity by Descent (IBD) or Principal Component Analysis (PCA).

+ __Select a cohort__: Choose a cohort from either the 1000 Genomes Project or the Human Genome Diversity Project.

+ __Select a resolution value__: Select one of the 50 resolution values. The resolution values span a logarithmic space from -2 to 2.

As a result of setting the previous options, you will obtain: 

+ __Resolution plot__: This plot displays the identified communities at various resolution values by genetic metric. At the bottom of the app you can activate the option "Similar colors" which enable communities of greater similarity (in a genetic sense) to be represented with more similar colors, promoting visual coherence. 

+ __Map__: In this map, cohorts are spatially linked to their respective sampling locations. Each pie chart corresponds to a cohort sourced from either the 1000 Genomes Project or The Human Genome Diversity Project.

+ __Pie chart__: The pie chart illustrates the distribution of individuals within a chosen cohort across various network communities.

+ __Community network__: This network is constructed by averaging the positions of individuals within each community in the network of individuals. The node sizes are proportional to the respective community sizes. Edges indicate the density of connections between these communities.

Note that colors in the resolution plot correspond to those displayed in the accompanying map.

If you want to know more about how the shiny app works and the interpretation of the results you can consult the [paper]() and the [tutorials]() made by the research team 

## Run the pipeline on the command-line

+ __Clone the repository__: First, you need to clone the GitHub repository. You can do this by using the __'git clone'__ command.
  <center>
    
  For example:
```
  git clone https://github.com/mariajpalma/GG-NC.git

 ``` 
  </center>
  
+ __Change to the newly cloned directory__: Once the repository is cloned, access the directory that was created during cloning using the __'cd'__ command.
  <center>
```    
  cd GG-NC
```  
  </center>

 + __Execute the '.sh' file__: Once you are in the repository directory, you can execute the __'.sh'__ file using the __'bash'__ command followed by the file name and the necessary parameters. You can see the required parameters in the ParametersV011223.sh file.
   <center>
    
  For example:
```
  bash ParametersV011223.sh -k <param1_value> -p <param2_value> -d <param3_value> -i <param4_value> -m <param5_value> -s <param6_value> -l <param7_value> -n <param8_value>
 ```
 Where each parameter refers to: 
 
 -k \<param1\>:   The genetic metric to consider. Kind: IBD, PCA or GRM
 
 -p \<param2\>:   Path of your files
 
 -d \<param3\>:   Name of your data files
 
 -i \<param4\>:    Name of your info file
 
 -m \<param5\>:   Maximum value
 
 -s \<param6\>:   Number of steps inside lambda threshold
 
 -l \<param7\>:   Lambda or resolution value to explore (logarithmic space from -2 to 2)
 
 -n \<param8\>:   Info file to generate shiny input
   </center>
   
 + Make sure to replace \<paramX_value\> with the specific values you want to pass as arguments to the script.
 
  For example:
```
  bash ParametersV011223.sh -k IBD -p /path/to/files -d data_file_name -i info_file_name -m 100 -s 5 -l 0.5 -n shiny_info_file
 ``` 
   </center>



With these steps, you should be able to run the shell script found in the GitHub repository from your system terminal. Remember that you will need to have Git installed on your system to clone the repository from GitHub.

## Authors 

María J. Palma, Centro de Ciencias Genómicas, UNAM, Cuernavaca, México.

Yuridia Posadas, Centro de Ciencias Genómicas, UNAM, Cuernavaca, México.

Claudia Quiroz, Escuela Nacional de Antropología e Historia, CDMX, México.

Brenda E. López, Centro de Ciencias Genómicas, UNAM, Cuernavaca, México. 

Anna Lewis, E J Safra Center for Ethics, Harvard University, Cambridge, MA.

Tina Lasisi, Univ. of Michigan, Ann Arbor, MI.

Kevin A. Bird, Univ. of California, Davis, Davis, CA.

Arslan Zaidi, Department of Genetics, cell biology, and Development, University of Minnesota, Twin Cities, MN/ Institute of Health Informatics, University of Minnesota, Twin Cities, MN.

Mashaal Sohail, Centro de Ciencias Genómicas, UNAM, Cuernavaca, México.

## Contact

mashaal@ccg.unam.mx
