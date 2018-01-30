-----------------------------------------------------------------------------------------------------------------------
------------  bcdLL1_alsls archive----------------------------------------------------------------------------------------
------------  Block-Component-Decomposition in rank-(L,L,1) terms of a third order tensor via ALS coupled with Line Search
------------  @Copyright Dimitri Nion  -------------------------------------------------------------------------------
------------- Released May 2010 ---------------------------------------------------------------------------------------
------------  Feedback dimitri.nion@gmail.com -------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------

Content of the pack: 
- Main function: bcdLL1_alsls.m (all subfunctions nedeed are strapped in bcdLL1_alsls.m)
- demo1, demo2, demo3, demo4: to show how to use the main function
- bcdLL1_init.m : this is the same as the bcdLL1_init.m subfunction strapped in bcdLL1_alsls.m
It has been copied in the folder because it is called by demo1 and demo3
- solve_blockperm.m : remove block-permutation and block-rotation ambiguity 
before computation of the error on the loading matrices estimates.