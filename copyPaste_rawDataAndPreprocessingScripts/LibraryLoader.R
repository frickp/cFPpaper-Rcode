# Function to load a package, or install it if it is not found.

# Example usage:	getLib("myLibrary"); getLib("myLibrary")
# Note, using twice consecutively allows the package to be simultaneously installed and loaded

getLib = function(myLib)
{
	if(eval(parse(text=paste0('library(',myLib,',logical.return=T)'))))
		{print(paste('package','found'))} else {	
		{
			install.packages(myLib,repos='http://cran.us.r-project.org')
			print('Package installed. Please re-run script to load library into memory.')
		}
	}
}
