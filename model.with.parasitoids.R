########## parameters/variables ####################

# f is the proportion of female offspring for host
# v is the daily per capita fecundity for host
# s is survivourship
# t is stage duration
# q is an exponential decay factor

# F is the daily per capita production of female offspring

# P is the probability of survival and remaining within the same life stage (P=s^1/t(1-1/t))
# G is the probability of survival and transitioning to the next life stage (P=s^1/t(1/t))

##Paropsis charybdis
#E1 are young eggs
#E2 are old eggs
#L are larvae
#P are pupae
#PreA are pre-sexually mature adults
#A are adults

##parasitiods
#EN1 are within-egg life stages Enoggera nassaui
#EN2 are adult Enoggera nassaui
#NE1 are within-egg life stages of Neopolycystus insectifurax
#NE2 are adult Neopolycystus insectifurax


### set parameters ####

#survivourship
E1s <- 0.9523
E2s <- 0.9523
Ls <- 0.45
Ps <- 0.715
PreAs <- 1

EN1s <- 0.8241208

NE1s <- 0.8241208


#duration
E1t <- 2
E2t <- 3
Lt <- 17.7
Pt <- 12.5
PreAt <- 15
At <- 60

EN1t <- 9
EN2t <- 16.7

NE1t <- 11
NE2t <- 16.7

#overwintering survival parameter for Enog adults
EN2.ovwint <- 0.4


#daily per capita production of female offspring paropsis
E1f <- 0.5
E1v <- 17.3


#mean #eggs laid per parasitoid
EnogB <- 9.187624551
NeopB <- 9.187624551


## Probability of surviving and remaining within the same life stage, P
E1P <- (E1s^(1/E1t))*(1-(1/E1t))
E2P <- (E2s^(1/E2t))*(1-(1/E2t))
LP <- (Ps^(1/Pt))*(1-(1/Pt))
PP <- (Ps^(1/Pt))*(1-(1/Pt))
PreAP <- (PreAs^(1/PreAt))*(1-(1/PreAt))
AP <- (1-(1/At))^2

EN1P <- (EN1s^(1/EN1t))*(1-(1/EN1t))
EN2P <- (1-(1/EN2t))^2

NE1P <- (NE1s^(1/NE1t))*(1-(1/NE1t))
NE2P <- (1-(1/NE2t))^2

## Probability of surviving and moving onto the next life stage, G
E1G <- (E1s^(1/E1t))*(1/E1t)
E2G <- (E2s^(1/E2t))*(1/E2t)
LG <- (Ls^(1/Lt))*(1/Lt)
PG <- (Ps^(1/Pt))*(1/Pt)
PreAG <- (PreAs^(1/PreAt))*(1/PreAt) 

EN1G <- (EN1s^(1/EN1t))*(1/EN1t)


NE1G <- (NE1s^(1/NE1t))*(1/NE1t)



# exponential decay factor
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

EN1PopSize <- numeric(TotalDays)
EN2PopSize <- numeric(TotalDays)

NE1PopSize <- numeric(TotalDays)
NE2PopSize <- numeric(TotalDays)

####### Set initial population sizes #######
E1PopNow <- 0
E2PopNow <- 0
LPopNow <- 0
PPopNow <- 1898.331
PreAPopNow <- 2627.761
APopNow <- 0

EN1PopNow <- 0
EN2PopNow <- 0

NE1PopNow <- 0
NE2PopNow <- 0


############ Simulation loop parasitoids only ############
for (day in 1:TotalDays) { 
  if (day == 60) {
    #emergence of e. nassaui
    EN1PopNow <- 5140.46
    EN2PopNow <- 0
  }
  if (day == 120) {
    #emergence of n. insectifurax
    NE1PopNow <- 6746.462
    NE2PopNow <- 0
  }
  if (day >= 100) {
    #implementing overwintering of e. nassaui
    EnogB <- 0
    EN1G <- 0
    EN1P <- 1
    EN2.ovwint <- 0.4
  } else {
    EnogB <- 9.187624551
    EN1G <- (EN1s^(1/EN1t))*(1/EN1t)
    EN1P <- (EN1s^(1/EN1t))*(1-(1/EN1t))
    EN2.ovwint <- 1
  }
  # effect of sprays
  E1COR <- E1PopNow
  E2COR <- E2PopNow
  LCOR <- LPopNow
  PCOR <- PPopNow
  PreACOR <- PreAPopNow
  ACOR <- APopNow
  
  EN1COR <- EN1PopNow
  EN2COR <- EN2PopNow
  
  NE1COR <- NE1PopNow
  NE2COR <- NE2PopNow
  
  # Calculate egg production with exponential decay
  Fecund <- E1f * E1v * exp(-q * ACOR) * ACOR
  
  # number of host eggs that escaped parasitism from enog/were parasitised
  E1W <-  (ifelse(E1COR==0, 1, exp(-EnogB*EN2COR/E1COR))) * E1COR
  EN1W <- (ifelse(E1COR==0, 1, (1-exp(-EnogB*EN2COR/E1COR)))) * E1COR
  
  # number of host eggs that escaped parasitism from neop/were parasitised
  E1X <- ifelse(E1W==0,1, exp(-NeopB*NE2COR/E1W)) * E1W
  NE1X <- ifelse(E1W==0,1, (1-exp(-NeopB*NE2COR/E1W))) * E1W
  
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
  
  ##Keeping track of enog that move on
  EN1Mature <- (EN1COR + EN1W) * EN1G 
  
  ##Keeping track of enog that stay
  EN1PopNext <- (EN1COR + EN1W) * (1 - EN1G) * EN1P 
  EN2PopNext <- (EN2COR + EN1Mature) * EN2P * EN2.ovwint
  
  ##Keeping track of neop that move on
  NE1Mature <- (NE1COR + NE1X) * NE1G 
  
  ##Keeping track of neop that stay
  NE1PopNext <- (NE1COR + NE1X) * (1 - NE1G) * NE1P
  NE2PopNext <- (NE2COR + NE1Mature) * NE2P
  
  # Save results to vectors
  E1PopSize[day] <- E1PopNext
  E2PopSize[day] <- E2PopNext
  LPopSize[day] <- LPopNext
  PPopSize[day] <- PPopNext
  PreAPopSize[day] <- PreAPopNext
  APopSize[day] <- APopNext
  
  EN1PopSize[day] <- EN1PopNext
  EN2PopSize[day] <- EN2PopNext
  
  NE1PopSize[day] <- NE1PopNext
  NE2PopSize[day] <- NE2PopNext
  
  # turn next into now
  E1PopNow <- E1PopNext
  E2PopNow <- E2PopNext
  LPopNow <- LPopNext
  PPopNow <- PPopNext
  PreAPopNow <- PreAPopNext
  APopNow <- APopNext
  
  EN1PopNow <- EN1PopNext
  EN2PopNow <- EN2PopNext
  
  NE1PopNow <- NE1PopNext
  NE2PopNow <- NE2PopNext
}



########## no parasitoids or sprays x ##########

for (day in 1:TotalDays) { 
  # changing popnow to "effect of sprays" (there are no sprays)
  E1COR <- E1PopNow
  E2COR <- E2PopNow
  LCOR <- LPopNow
  PCOR <- PPopNow
  PreACOR <- PreAPopNow
  ACOR <- APopNow
  
  EN1COR <- EN1PopNow
  EN2COR <- EN2PopNow
  
  NE1COR <- NE1PopNow
  NE2COR <- NE2PopNow
  
  # Calculate egg production with exponential decay
  Fecund <- E1f * E1v * exp(-q * ACOR) * ACOR
  
  
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
  
  # Save results to vectors
  E1PopSize[day] <- E1PopNext
  E2PopSize[day] <- E2PopNext
  LPopSize[day] <- LPopNext
  PPopSize[day] <- PPopNext
  PreAPopSize[day] <- PreAPopNext
  APopSize[day] <- APopNext
  
  # turn next into now
  E1PopNow <- E1PopNext
  E2PopNow <- E2PopNext
  LPopNow <- LPopNext
  PPopNow <- PPopNext
  PreAPopNow <- PreAPopNext
  APopNow <- APopNext
}
