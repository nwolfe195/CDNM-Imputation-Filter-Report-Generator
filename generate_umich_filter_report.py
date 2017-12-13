import os
from jinja2 import Environment, FileSystemLoader
from optparse import OptionParser
import hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import scoreatpercentile
import shutil

PATH = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_ENVIRONMENT = Environment(
        autoescape=False,
        loader=FileSystemLoader(os.path.join(PATH, 'templates')),
        trim_blocks=False)

def main():
        create_index_html()

def create_index_html():
	# Set the template file to be used.  Jinja2 will be default look inside the templates directory
        temp = "generate_umich_filter_report_template.html"
	# Read in command line arguments.  Example format -f /proj/regeps/regep00/studies/VDAART/data/gwas/2016_VDAART_Mother_GWAS
        parser = OptionParser()
        parser.add_option("-f", "--folder", action="store", dest="folder")
        parser.add_option("-n", "--name", action="store", dest="name")
        options,args = parser.parse_args()
	# Read directory path and file name into variables
        if not options.folder.endswith('/'):
                options.folder = "%s%s" % (options.folder, '/')
        if not options.name:
                options.name = options.folder.split('/')[-2]

	# Set output directory
	out_dir = "%s_output" % (options.name)
	
	# Delete existing output directory and reinitialize if necessary
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)
	os.makedirs(out_dir)        

	# Output filename, within the output directory
	fname = "%s/%s_summary.html" % (out_dir, options.name)

	# Files of interest for the report
        files = ["Chromosome", "Exclude", "Force-Allele1", "ID", "Position", "Strand-Flip"]

	# Various methods to call and get relevent data
        percentage = get_percentage(len(files), options.folder, options.name, files)
        md5sums_files = get_md5sums_files(len(files), options.name,files)
        log_lines = get_log_lines(options.folder, options.name)
        md5sums = get_md5sums(len(files), options.folder, md5sums_files)
        FreqPlots, Outliers = get_FreqPlots(options.folder, options.name, out_dir)

	# Call render_template, passing all of the method created data
        with open(fname, 'w') as f:
                html = render_template(temp, files, options.name, percentage, options.folder, log_lines, md5sums_files, md5sums, FreqPlots, Outliers)
                f.write(html)

def get_percentage(file_num, folder, name, files):
	# Initialize 2d list for output
        percentage = [[0]*(file_num+1) for _ in range(22)]
        # Iterate through each chromosome
	for i in range(1,23):
		# Get file path for preimputation file
                totalFile = "%schr%d/%s_chr%d.bim" % (folder, i, name, i)
		# Get total number of snps (1 per line)
                total = len([line.rstrip('\n') for line in open(totalFile)])
		# Set first column to the total from above
                percentage[i-1][0] = total
		# Iterate through each file type of interest
                for f in files:
			# Get filepath for individual file type
                        file = "%spre_imputation/chr%d/%s-%s_chr%d-HRC.txt" % (folder, i, f, name, i)
                        # Get number of snps (1 per line)
			number = float(len([line.rstrip('\n') for line in open(file)]))
                 	# Calculate percentage
			percent = (number/total)*100
			# Round to make it pretty
                        percent = "%.2f" % round(percent, 2)
			# Pass total number and percentage to the output 2D list
                        percentage[i-1][files.index(f)+1] = "%d(%s%%)" % (number, percent)
	# Pass the 2D list back to the main program
        return percentage

def get_md5sums_files(file_num, name, files):
	# Intialize 2D list to hold files that need md5 sums
        md5sums_files = [[0]*(file_num+2) for _ in range(22)]
        # Iterate through each chromosome
	for i in range(1,23):
                # Get pre imputation loaded in
		md5sums_files[i-1][0] = "chr%d/%s_chr%d.bim" % (i, name, i)
		# Iterate through files of interest
                for f in files:
			# Add filename of each file of interest
                        md5sums_files[i-1][files.index(f)+1] = "pre_imputation/chr%d/%s-%s_chr%d-HRC.txt" % (i, f, name, i)
		# Add final file of interest not on the file list (FreqPlot)
                md5sums_files[i-1][len(files)+1] = "pre_imputation/chr%d/%s-%s_chr%d-HRC.txt" % (i, "FreqPlot", name, i)
	# Return files to md5 sum
        return md5sums_files

def get_log_lines(folder, name):
	# Get example log file.  Same output for all of them, so arbitrarily choosing chromosome 1 to use
        log_file = "%spre_imputation/chr%d/%s-%s_chr%d-HRC.txt" % (folder, 1, "LOG", name, 1)
	# Read in the file
        with open(log_file) as f:
                log_lines = f.read().splitlines()
        return log_lines

def get_md5sums(file_num, folder, md5sums_files):
	# Initialize md5sum 2D list
        md5sums = [[0]*(file_num+2) for _ in range(22)]
	# Iterate through the 2D list
        for i in range(0,22):
                for j in range(0, (file_num+2)):
			# Construct current file
                        md5_file = "%s%s" % (folder, md5sums_files[i][j])
			# Read file in
                        with open(md5_file) as f:
                                file_content = f.read()
			# Hash the contents of the file into md5 sum
                        md5sums[i][j] = hashlib.md5(file_content).hexdigest()
	# Return md5sum 2D list
        return md5sums

def get_FreqPlots(folder, name, out_dir):
	# Initialize Freqplot and Outlier lists
        FreqPlots = [0]*22
        Outliers = [0]*22
	# Iterate through each chromosome
        for i in range(1,23):
		# Initialize address of frequency plots
                FreqPlots[i-1] = "%s/plot%d.png" % (out_dir, i)
		# Set file to get frequency data
                freq_file = "%spre_imputation/chr%d/%s-%s_chr%d-HRC.txt" % (folder, i, "FreqPlot", name, i)
		# Read in freq_file data
		with open(freq_file) as f:
                        data = [line.split() for line in f]
                # Make plots and output
		x = map(float, [item[1] for item in data])
                y = map(float, [item[2] for item in data])
                curr_plot = plt.figure()
                plt.plot(x,y,'.',figure=curr_plot)
                curr_plot.suptitle("Chromosome %d" % (i))
                curr_plot.savefig(FreqPlots[i-1])
                plt.close(curr_plot)
	
		# Calculate outliers, grab the top 10
                q1, q3 = np.percentile(x, [25,75])
                upper = q3+1.5*(q3-q1)
                outliers = [[x,el] for [x,el] in enumerate(x) if el>upper]
                outliers = sorted(outliers, key=lambda x: float(x[1]), reverse=True)[:10]
                for j in range(0,10):
                        outliers[j][0] = data[j][0]
                Outliers[i-1] = outliers
        return FreqPlots, Outliers

def render_template(template_filename, files, name, percentage, path, log_lines, md5sums_files, md5sums, freqplots, outliers):
        # List of variables to be passed
	templateVars = { "files" : files,
                         "dataset" : name,
                         "tabledata" : percentage,
                         "path" : path,
                         "log" : log_lines,
                         "md5" : md5sums_files,
                         "md5sums" : md5sums,
                         "plot_paths" : freqplots,
                         "outliers" : outliers
                        }
	# Render and return completed HTML file
        return TEMPLATE_ENVIRONMENT.get_template(template_filename).render(templateVars)

if __name__ == "__main__":
        main()

