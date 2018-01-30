----------------------------------------------------------------------------------------------------------------------
------------  bcdLM_alsls archive----------------------------------------------------------------------------------------
------------  Block-Component-Decomposition in rank-(Lr,Mr,.) terms of a third order tensor via ALS coupled with Line Search
------------  @Copyright Dimitri Nion  -------------------------------------------------------------------------------
------------- Released May 2010 ---------------------------------------------------------------------------------------
------------  Feedback dimitri.nion@gmail.com -------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------

Content of the pack: 
- Main function: bcdrLMr_alsls.m (all subfunctions nedeed are strapped in bcdLrMr_alsls.m)
- bcd_LrMr_init: to generate a random initialization
- demo: to show how to use the main function
- solve_blockperm.m : remove block-permutation and block-rotation ambiguity 
before computation of the error on the loading matrices estimates.