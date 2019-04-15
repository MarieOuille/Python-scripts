# Fringeezz

#Before starting, you will need to have a hdf5 file recorded with the Fringeezz (data_logger mode). 
#You will also need to know the name of the file (usually phase_log_date.hdf5) and the path of the file (usually in SalleNoire2Network).

#The file FRINGEEZZ.py then just needs to be executed with Python in a command prompt (I personnally use Console2 : python FRINGEEZZ.py)
#It will ask you where the file is located and what is the name of the file : just copy/paste the path and name in the command prompt

#OUTPUT :
#It automatically generates and saves 2 plots : 1) CEP as a function of time + histogram + standard deviation value and 2) power spectrum density + integrated phase noise.
#It also creates a basic txt file with 2 columns : CEP values (in rad, centered around 0) ; time (in min, starting at 0, with 1ms btw each value) 

#The 4 files starting with phase_log in this repository are just examples (The hdf5 file was recorded with the Fringeezz ; the png and txt files are then generated when executing the FRINGEEZZ.py file)

#The file screenshot_HDFview.png shows how data are recorded in the hdf5 format. You will find more explanations in the comments in the FINGEEZZ.py file 
