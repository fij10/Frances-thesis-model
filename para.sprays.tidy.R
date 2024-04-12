########## parameters/variables ####################

# f is the proportion of female offspring for host
# v is the daily per capita fecundity for host
# s is survivourship
# t is stage duration
# q is an expoNeopopopntial decay factor

# F is the daily per capita production of female offspring

# P is the probability of survival and remaining within the same life stage (P=s^1/t(1-1/t))
# G is the probability of survival and transitioning to the Next life stage (P=s^1/t(1/t))

##Paropsis charybdis
#E1 are young eggs
#E2 are old eggs
#L are larvae
#P are pupae
#PreA are pre-sexually mature adults
#A are adults

##parasitiods
#Enog1 are within-egg life stages Enoggera nassaui
#Enog2 are adult Enoggera nassaui
#Neop1 are within-egg life stages of Neopolycystus insectifurax
#Neop2 are adult Neopolycystus insectifurax


### set parameters ####

#survivourship
E1s <- 0.9523
E2s <- 0.9523
Ls <- 0.45
Ps <- 0.715
PreAs <- 1

Enog1s <- 0.8241208

Neop1s <- 0.8241208


#duration
E1t <- 2
E2t <- 3
Lt <- 17.7
Pt <- 12.5
PreAt <- 15
At <- 60

Enog1t <- 9
Enog2t <- 16.7

Neop1t <- 11
Neop2t <- 16.7

#overwintering survival parameter for Enog adults
Enog2.ovwint <- 0.4


#daily per capita production of female offspring paropsis
E1f <- 0.5
E1v <- 17.3


#mean #eggs laid per parasitoid
EnogogB <- 9.187624551
NeopopB <- 9.187624551


## Probability of surviving and remaining within the same life stage, P
E1P <- (E1s^(1/E1t))*(1-(1/E1t))
E2P <- (E2s^(1/E2t))*(1-(1/E2t))
LP <- (Ps^(1/Pt))*(1-(1/Pt))
PP <- (Ps^(1/Pt))*(1-(1/Pt))
PreAP <- (PreAs^(1/PreAt))*(1-(1/PreAt))
AP <- (1-(1/At))^2

Enog1P <- (Enog1s^(1/Enog1t))*(1-(1/Enog1t))
Enog2P <- (1-(1/Enog2t))^2

Neop1P <- (Neop1s^(1/Neop1t))*(1-(1/Neop1t))
Neop2P <- (1-(1/Neop2t))^2

## Probability of surviving and moving onto the Neopxt life stage, G
E1G <- (E1s^(1/E1t))*(1/E1t)
E2G <- (E2s^(1/E2t))*(1/E2t)
LG <- (Ls^(1/Lt))*(1/Lt)
PG <- (Ps^(1/Pt))*(1/Pt)
PreAG <- (PreAs^(1/PreAt))*(1/PreAt) 

Enog1G <- (Enog1s^(1/Enog1t))*(1/Enog1t)


Neop1G <- (Neop1s^(1/Neop1t))*(1/Neop1t)



# expoNeopntial decay factor
q <- 0.0005

# Set total number of days
TotalDays <- 183

# Initialize vectors to hold results
E1PopSize <- numeric(TotalDays)
E2PopSize <- numeric(TotalDays)
LPopSize <- numeric(TotalDays)
PPopSize <- numeric(TotalDays)
PreAPopSize <- numeric(TotalDays)
APopSize <- numeric(TotalDays)

Enog1PopSize <- numeric(TotalDays)
Enog2PopSize <- numeric(TotalDays)

Neop1PopSize <- numeric(TotalDays)
Neop2PopSize <- numeric(TotalDays)

######################## simulations #################
# spray survivorship

#survivorship due to sprays
#Spray oNeop
E1Css <- 0.95
E2Css <- 0.95
LCss <- 0
PCss <- 1
PreACss <- 0.77
ACss <- 0.77

Enog1Css <- 0.77
Enog2Css <- 0.5375

Neop1Css <- 0.77
Neop2Css <- 0.5375

#Spray two
E1Css <- 0.84
E2Css <- 0.84
LCss <- 0
PCss <- 1
PreACss <- 0.325
ACss <- 0.325

Enog1Css <- 1
Enog2Css <- 0

Neop1Css <- 1
Neop2Css <- 0

#Spray three
E1Css <- 0.00001
E2Css <- 0.00001
LCss <- 0.00001
PCss <- 1
PreACss <- 0.00001
ACss <- 0.00001

Enog1Css <- 0.00001
Enog2Css <- 0.00001

Neop1Css <- 0.00001
Neop2Css <- 0.00001

#application day

spray.day <- 30

#initial population sizes before overwintering + next simulation
PreAPopNow <-  1067.462
APopNow <-  1528.882
Enog1PopNow <-  4894.255
Neop1PopNow <- 6717.555


# overwintering survival probabilities
# adult p.c
ovwinter.pc <- (PreAPopNow + APopNow) * 0.5
# within-egg e.n
ovwinter.enog <- Enog1PopNow * 0.78
# within-egg n.i
ovwinter.neop <- Neop1PopNow * 0.78

# Set initial population sizes for sim
E1PopNow <- 0
E2PopNow <- 0
LPopNow <- 0
PPopNow <- 0
PreAPopNow <- ovwinter.pc
APopNow <- 0

Enog1PopNow <- 0
Enog2PopNow <- 0

Neop1PopNow <- 0
Neop2PopNow <- 0


##### loop #####

for (day in 1:TotalDays) { 
  if (day == 60) {
    #emergence of e. nassaui
    Enog1PopNow <- ovwinter.enog
    Enog2PopNow <- 0
  }
  if (day == 120) {
    #emergence of n. insectifurax
    Neop1PopNow <- ovwinter.neop
    Neop2PopNow <- 0
  }
  if (day >= 100) {
    #implemEnogting overwintering of e. nassaui
    EnogogB <- 0
    Enog1G <- 0
    Enog1P <- 1
    Enog2.ovwint <- 0.4
  } else {
    EnogogB <- 9.187624551
    Enog1G <- (Enog1s^(1/Enog1t))*(1/Enog1t)
    Enog1P <- (Enog1s^(1/Enog1t))*(1-(1/Enog1t))
    Enog2.ovwint <- 1
  }
  if (day == spray.day) {
    #survivourship due to sprays
    E1Css <- 0.84
    E2Css <- 0.84
    LCss <- 0
    PCss <- 1
    PreACss <- 0.325
    ACss <- 0.325
    
    Enog1Css <- 1
    Enog2Css <- 0
    
    Neop1Css <- 1
    Neop2Css <- 0
  } else {
    E1Css <- 1
    E2Css <- 1
    LCss <- 1
    PCss <- 1
    PreACss <- 1
    ACss <- 1
    
    Enog1Css <- 1
    Enog2Css <- 1
    
    Neop1Css <- 1
    Neop2Css <- 1
  }
  # effect of sprays
  E1COR <- E1Css * E1PopNow
  E2COR <- E2Css * E2PopNow
  LCOR <- LCss * LPopNow
  PCOR <- PCss * PPopNow
  PreACOR <- PreACss * PreAPopNow
  ACOR <- ACss * APopNow
  
  Enog1COR <- Enog1Css * Enog1PopNow
  Enog2COR <- Enog2Css * Enog2PopNow
  
  Neop1COR <- Neop1Css * Neop1PopNow
  Neop2COR <- Neop2Css * Neop2PopNow
  
  # Calculate egg production with exponential decay
  Fecund <- E1f * E1v * exp(-q * ACOR) * ACOR
  
  # number of host eggs that escaped parasitism from Enog/were parasitised
  E1W <-  (ifelse(E1COR==0, 1, exp(-EnogogB*Enog2COR/E1COR))) * E1COR
  Enog1W <- (ifelse(E1COR==0, 1, (1-exp(-EnogogB*Enog2COR/E1COR)))) * E1COR
  
  # number of host eggs that escaped parasitism from Neop/were parasitised
  E1X <- ifelse(E1W==0,1, exp(-NeopopB*Neop2COR/E1W)) * E1W
  Neop1X <- ifelse(E1W==0,1, (1-exp(-NeopopB*Neop2COR/E1W))) * E1W
  
  # Update population sizes
  ##this is keeping track of those that move to the next life stage
  E1Mature <- (Fecund + E1X) * E1G
  E2Mature <- (E2COR) * E2G
  LMature <- (LCOR) * LG
  PMature <- (PCOR) * PG
  PreAMature <- (PreACOR) * PreAG
  
  ##this is keeping track of those that are staying in the same life stage 
  E1PopNext <- (Fecund + E1X) * (1 - E1G) * E1P
  E2PopNext <- (E2COR * (1 - E2G) + E1Mature) * E2P
  LPopNext <- (LCOR * (1 - LG) + E2Mature) * LP
  PPopNext <- (PCOR * (1 - PG) + LMature) * PP
  PreAPopNext <- (PreACOR * (1-PreAG) + PMature) * PreAP
  APopNext <- (ACOR * (1 - AG) + PreAMature) * AP
  
  ##Keeping track of Enogog that move on
  Enog1Mature <- (Enog1COR + Enog1W) * Enog1G 
  
  ##Keeping track of Enogog that stay
  Enog1PopNext <- (Enog1COR + Enog1W) * (1 - Enog1G) * Enog1P 
  Enog2PopNext <- (Enog2COR + Enog1Mature) * Enog2P * Enog2.ovwint
  
  ##Keeping track of Neopop that move on
  Neop1Mature <- (Neop1COR + Neop1X) * Neop1G 
  
  ##Keeping track of Neopop that stay
  Neop1PopNext <- (Neop1COR + Neop1X) * (1 - Neop1G) * Neop1P
  Neop2PopNext <- (Neop2COR + Neop1Mature) * Neop2P
  
  # Save results to vectors
  E1PopSize[day] <- E1PopNext
  E2PopSize[day] <- E2PopNext
  LPopSize[day] <- LPopNext
  PPopSize[day] <- PPopNext
  PreAPopSize[day] <- PreAPopNext
  APopSize[day] <- APopNext
  
  Enog1PopSize[day] <- Enog1PopNext
  Enog2PopSize[day] <- Enog2PopNext
  
  Neop1PopSize[day] <- Neop1PopNext
  Neop2PopSize[day] <- Neop2PopNext
  
  # turn Neopxt into now
  E1PopNow <- E1PopNext
  E2PopNow <- E2PopNext
  LPopNow <- LPopNext
  PPopNow <- PPopNext
  PreAPopNow <- PreAPopNext
  APopNow <- APopNext
  
  Enog1PopNow <- Enog1PopNext
  Enog2PopNow <- Enog2PopNext
  
  Neop1PopNow <- Neop1PopNext
  Neop2PopNow <- Neop2PopNext
}


