README
This python script produces a report detailing pre-imputation cleaning.
Run from the command line, the script accepts up to 2, but at least 1 command
line argument.  The required one, -f or --file, is a full filepath to the 
location containing chr directories.  The 2016_VDAART_Mother_GWAS data set 
for example has the path /proj/regeps/regep00/studies/VDAART/data/gwas/2016_VDAART_Mother_GWAS/
The second argument, -n or --name, is the name of the data set.  The 2016_VDAART_Mother_GWAS
data set has the name 2016_VDAART_Mother_GWAS.  This  argument is optional, as
the script will extract the name from the filepath.  Only use if the desired 
data set name is different from its filepath.

Using the Jinja2 template
This script uses the Jinja2 package to construct the html output page.  The template file 
generate_umich_filter_report_template.html must be placed in a directory named templates, 
stored in the same directory the script is being run from.

Output
The script will output an HTML file, named <name>_summary.html.  The 2016_VDAART_Mother_GWAS
data set would be named 2016_VDAART_Mother_GWAS_output.html.  In addition, one .png file 
will be created per chromosome, which the html file will need to properly display the 
plots contained within it.
Both the HTML file and the graphs will all be found in a directory named after the data
set used to construct the report, ie 2016_VDAART_Mother_GWAS_output.
