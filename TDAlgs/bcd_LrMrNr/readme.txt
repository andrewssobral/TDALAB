-----------------------------------------------------------------------------------------------------------------------
------------  bcdLMN_alsls archive----------------------------------------------------------------------------------------
------------  Block-Component-Decomposition in rank-(Lr,Mr,Nr) terms of a third order tensor via ALS coupled with Line Search
------------  @Copyright Dimitri Nion  -------------------------------------------------------------------------------
------------- Released May 2010 ---------------------------------------------------------------------------------------
------------  Feedback dimitri.nion@gmail.com -------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------

Content of the pack: 
- Main function: bcdLrMrNr_alsls.m (all subfunctions nedeed are strapped in bcdLrMrNr_alsls.m)
- demo1, demo2: to show how to use the main function
- tmprod: tensor-matrix product, called by the demo file to generate the tensor to be decomposed.
- bcdLrMrNr_init.m : this is the same as the bcdLrMrNr_init.m subfunction strapped in bcdLrMrNr_alsls.m
It has been copied in the folder because it is called by the demo script.
- solve_blockperm.m : remove block-permutation and block-rotation ambiguity 
before computation of the error on the loading matrices estimates.