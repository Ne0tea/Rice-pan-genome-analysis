Script for pan TE database building
  
##00_command_save  
terimal command used for pipeline  
  
##01_repeatmodeler.pbs  
batch repeatmodeler submit script
  
##02_ltr_class  
02_Ltr_and_classfier.pl used in repeatmodeler request for different memroy and thread, so seprate script was used for TE classifying.  
  
##03_deepTE_classfy.pbs  
For unclassfied TE, further classified in DeepTE  

##04_repeatmask.pbs  
batch repeatmask submit script
