netassess-new
=============

This is Mike Rizzo's netassess package, updated to work in new versions of R. It is still a work in progress. 

**NOTE**

I have removed the requirement for ```gpclib``` from the package. After doing some research, this appears to be an old package
with a non-open license, it was a dependency for another package used by ```netassess``` and has since been replaced with a package called ```rgeos```. I have updated the dependencies of ```netassess-new``` to reflect this and hopefully this fixes the problem. If it doesn't, Nathan is probably correct that you will need to install Rtools in order to compile the source package version of ```gpclib```. You can install Rtools from this link:

http://cran.r-project.org/bin/windows/Rtools/

If you discover that you still need ```gpclib```, once you have Rtools installed you can install it as before:

```install.packages("gpclib", type = "source")```

You can then install netassess-new from github with the following command:

```devtools::install_github("netassess-new", "ebailey78")```

If these instructions don't work you can contact me at ebailey@idem.in.gov.
