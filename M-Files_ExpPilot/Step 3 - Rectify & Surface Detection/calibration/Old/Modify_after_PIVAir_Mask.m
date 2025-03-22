function [XPIV_LFV_Surface,PIV_LFV_Surface,PIV_Surface] = Modify_after_PIVAir_Mask(XPIV_LFV_Surface,PIV_LFV_Surface,PIV_Surface)

XPIV_LFV_Surface = XPIV_LFV_Surface-40;
PIV_LFV_Surface = PIV_LFV_Surface-10;
PIV_Surface = PIV_Surface(41:4066)-10;