# Check whether the two functions produce similar estimates

## Disease setting
a <- intEffectAlpha(N = 1e5,
         setting = "Disease",
         eta = rep(0.1, 3),
         nu = rep(1.1, 3),
         beta_L0_D = 1,
         beta_L0_L = 1,
         beta_L_D = 1,
         beta_A0_D = 0,
         beta_A0_L = 0,
         beta_L0_C = 0,
         beta_A0_C = 0,
         beta_L_C = 0,
         cens = 0,
         years_lost = FALSE)

b <- intEffectAlphaDisease(N = 10^5,
                     beta_L0_D = 1,
                     beta_L0_L = 1,
                     beta_L_D = 1,
                     beta_A0_D = 0,
                     beta_A0_L = 0,
                     years_lost = FALSE)


## Drop In setting
intEffectAlpha(N = 1000, alpha = 0.5, tau = 5, years_lost = FALSE, a0 = 1, setting = "Drop In",
         beta_L_A = 1,
         beta_L_Z = 2,
         beta_L_D = 1.5,
         beta_L_C = 0,
         beta_A_L = -0.5,
         beta_A_Z = -0.5,
         beta_A_D = -1,
         beta_A_C = 0,
         beta_Z_L = -1,
         beta_Z_A = 0,
         beta_Z_D = -1,
         beta_Z_C = 0,
         beta_L0_L = 1,
         beta_L0_A = 1,
         beta_L0_Z = 1,
         beta_L0_D = 1,
         beta_L0_C = 0,
         beta_A0_L = -1.5,
         beta_A0_A = 0,
         beta_A0_Z = 0,
         beta_A0_D = -2,
         beta_A0_C = 0,
         cens = 0)

intEffectAlphaDropIn(N = 10^5,  beta_L_A = 1,
               beta_L_Z = 2,
               beta_L_D = 1.5,
               beta_L_C = 0,
               beta_A_L = -0.5,
               beta_A_Z = -0.5,
               beta_A_D = -1,
               beta_A_C = 0,
               beta_Z_L = -1,
               beta_Z_A = 0,
               beta_Z_D = -1,
               beta_Z_C = 0,
               beta_L0_L = 1,
               beta_L0_A = 1,
               beta_L0_Z = 1,
               beta_L0_D = 1,
               beta_L0_C = 0,
               beta_A0_L = -1.5,
               beta_A0_A = 0,
               beta_A0_Z = 0,
               beta_A0_D = -2,
               beta_A0_C = 0,
               lower = 10^(-400),
               upper = 10^3)

## Treatment setting
intEffectAlphaTreat(N = 10^5)
intEffectAlpha(N = 10^5, setting = "Treatment", cens = 0)

