-----------------------------------------------------------------------------------------------------------------------
------------  bcdLMN_alsls archive----------------------------------------------------------------------------------------
------------  Block-Component-Decomposition in rank-(L,M,N) terms of a third order tensor via ALS coupled with Line Search
------------  @Copyright Dimitri Nion  -------------------------------------------------------------------------------
------------- Released May 2010 ---------------------------------------------------------------------------------------
------------  Feedback dimitri.nion@gmail.com -------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------

Content of the pack: 
- Main function: bcdLMN_alsls.m (all subfunctions nedeed are strapped in bcdLMN_alsls.m)
- demo1, demo2: to show how to use the main function
- tmprod: tensor-matrix product, called by the demo file to generate the tensor to be decomposed.
- bcdLMN_init.m : this is the same as the bcdLMN_init.m subfunction strapped in bcdLMN_alsls.m
It has been copied in the folder because it is called by the demo script.
- solve_blockperm.m : remove block-permutation and block-rotation ambiguity 
before computation of the error on the loading matrices estimates.