# Integrated analysis of spotted seal fish biomass consumption
# Paul Conn, Last updated 1/22/26

library(raster) 
library(sf)
library(ggplot2)
set.seed(123456)
#Abundance of seals

N_dist = exp(rnorm(10000,12.45542,sqrt(0.00556712)))  #total abundance of spotted seals (from calculate_seal_abundance.R)

#Spatial distribution

tif_winter=raster("./Data/Combined_spotted_winter_cell_prob.tif")
tif_summer=raster("./Data/Combined_spotted_summer_cell_prob.tif")
winter_UD_poly=rasterToPolygons(tif_winter)
summer_UD_poly=rasterToPolygons(tif_summer)
winter_UD = st_centroid(st_as_sf(winter_UD_poly))
summer_UD = st_centroid(st_as_sf(summer_UD_poly))
winter_UD = st_transform(winter_UD,4326)
summer_UD = st_transform(summer_UD,4326)
lat_winter = st_coordinates(winter_UD)[,"Y"]  
lat_summer = st_coordinates(summer_UD)[,"Y"]
Strait_lat = 65.8  #this corresponds to little diomede

Summer_prop_chukchi = sum(summer_UD$Combined_spotted_summer_cell_prob[which(lat_summer>Strait_lat)]) 
Winter_prop_chukchi = sum(winter_UD$Combined_spotted_winter_cell_prob[which(lat_winter>Strait_lat)]) 


#Age-sex composition
load("./Data/Stable_stage_matrix.RData")
Age_proportions = Stable[4,]  #from Conn and Trukhanova 2022; 4th row is for spotted seals
SA_proportions = 0.5*matrix(Age_proportions,40,2)  #sex on column, age on row, 50/50 sex ratio

#bar plot of age-sex composition
plot_df = data.frame("Age"=c(0:39),"Proportion"=Age_proportions)
png("Age_proportions.png",width = 8, height = 8, units='in', res = 600)
  ggplot(plot_df)+geom_bar(aes(x=Age,y=Proportion),stat="identity")+
    theme(text = element_text(size=14))+xlab("Age (years)")+ylab("Proportion of spotted seals")
dev.off()

#Number of seals by age, sex, region, season, posterior replicate
#not currently including inter-year mortality
N_asrtk = array(NA,dim=c(40,2,2,2,10000))
for(ireg in 1:2){
  for(iseas in 1:2){
    for(irep in 1:10000){
      N_asrtk[,,ireg,iseas,irep]=SA_proportions*N_dist[irep]
    }
  }
}
N_asrtk[,,1,1,] = (1-Summer_prop_chukchi)*N_asrtk[,,1,1,]  #Bering, summer
N_asrtk[,,1,2,] = (1-Winter_prop_chukchi)*N_asrtk[,,1,2,]  #Bering, winter
N_asrtk[,,2,1,] = Summer_prop_chukchi*N_asrtk[,,2,1,]  #Chukchi, summer
N_asrtk[,,2,2,] = Winter_prop_chukchi*N_asrtk[,,2,2,]  #Chukchi, winter

#table for numbers at age
N_mean = apply(N_asrtk,c(1:4),'mean')
N_mean = apply(N_mean,c(1,3,4),'sum') #sum over sex
for(ir in 1:2){
  for(it in 1:2){
    N_mean[2,ir,it]=sum(N_mean[2:40,ir,it])
  }
}
N_mean=N_mean[1:2,,]
N_df = data.frame("Region_Season"=c("Bering/Summer","Chukchi/Summer",
                                    "Bering/Winter","Chukchi/Winter"),
                  "Pup" = as.vector(N_mean[1,,]),
                  "Non-pup"=as.vector(N_mean[2,,]))
N_df$Total = N_df[,2]+N_df[,3]

write.csv(N_df,file="N_age_region.csv")


####  Total energy consumption  ####

load('./Data/Seal_wt_age.RData')  #produced with seal_mass_at_age.Rmd
load('./Data/Energy_predictions.RData') #produced with captiveFoodIntake_massSpecific_kcal-kg.Rmd

#we'll calculate energy consumption per age class and month to start with
#to allow for different masses of seals within a season (important for young seals)
#before summing to get seasonal requirements 
#note that in the energy intake data, juvenile = ages 2,3,4; adult = 5+
#note also that actual month = calendar month - 1 in "Energy_predictions"
#but that in our integrated analysis the first month is May (first month of open water period) 

Month_order = c(5:12,1:4)
Energy_kcal = array(0,dim=c(40,12,10000))  #age, month (May = 1), sex, posterior replicate
#pups
Energy_kcal[1,2:12,] = 110  #this is the kcal/kg value assumed for pups in the open water period; 
#note no calories in first month (this is accounted for via adult female lactation)

#yearlings
Energy_kcal[2,1:12,] = Energy_predictions$pred_matrix[13:24,][Month_order]
#juveniles
for(iage in 2:4){
  Energy_kcal[iage,1:12,] = Energy_predictions$pred_matrix[25:36,][Month_order]
}
#adults
for(iage in 5:40){
  Energy_kcal[iage,1:12,] = Energy_predictions$pred_matrix[1:12,][Month_order]
}

#multiply by days in month
month30 = c(2,5,7,12)  #first month is may
month31 = c(1,3,4,6,8,9,11)
for(imonth in 1:length(month30)){
  Energy_kcal[,month30[imonth],]=Energy_kcal[,month30[imonth],]*30
}
for(imonth in 1:length(month31)){
  Energy_kcal[,month31[imonth],]=Energy_kcal[,month31[imonth],]*31
}
Energy_kcal[,10,]=Energy_kcal[,10,]*28  #Feb.

Energy_kcal_m = Energy_kcal_f = Energy_kcal #sex specific


#mutiply by mass-at-age
counter=1
for(iage in 1:40){
  for(imonth in 1:12){
    Energy_kcal_m[iage,imonth,] = Energy_kcal_m[iage,imonth,]*Wt_age$male[counter,]
    Energy_kcal_f[iage,imonth,] = Energy_kcal_f[iage,imonth,]*Wt_age$fem[counter,]
    counter=counter+1
  }
}

#cost of lactation
Energy_kcal_f[4,12,]=Energy_kcal_f[4,12,] + 0.333*382320  #1/3 of 3 year old females sexually mature
Energy_kcal_f[5:40,12,]=Energy_kcal_f[5:40,12,]+ 382320


#Now we posterior predictions for energy requirements by age (year, month) per seal. Next
#we need to aggregate these into age and season, and multiply by the number of seals in each
#age class

Kcal_tot_asrtk = array(NA,dim=c(40,2,2,2,10000))

for(iage in 1:40){
  for(ireg in 1:2){
    #open water
    Kcal_tot_asrtk[iage,1,ireg,1,] = Energy_kcal_f[iage,1,]
    Kcal_tot_asrtk[iage,2,ireg,1,] = Energy_kcal_m[iage,1,]
    for(imo in 2:7){  #sum over open water month
      Kcal_tot_asrtk[iage,1,ireg,1,] = Kcal_tot_asrtk[iage,1,ireg,1,]+Energy_kcal_f[iage,imo,]
      Kcal_tot_asrtk[iage,2,ireg,1,] = Kcal_tot_asrtk[iage,2,ireg,1,]+Energy_kcal_m[iage,imo,]
    }
    #winter
    Kcal_tot_asrtk[iage,1,ireg,2,] = Energy_kcal_f[iage,8,]
    Kcal_tot_asrtk[iage,2,ireg,2,] = Energy_kcal_m[iage,8,]
    for(imo in 9:12){  #sum over winter month
      Kcal_tot_asrtk[iage,1,ireg,2,] = Kcal_tot_asrtk[iage,1,ireg,2,]+Energy_kcal_f[iage,imo,]
      Kcal_tot_asrtk[iage,2,ireg,2,] = Kcal_tot_asrtk[iage,2,ireg,2,]+Energy_kcal_m[iage,imo,]
    }
    
  }
}
Kcal_tot_asrtk = Kcal_tot_asrtk * N_asrtk

## Now collapse into two age groups (pups / non-pups), combine sexes
Kcal_tot_artk = array(NA,dim=c(2,2,2,10000))
for(ir in 1:2){
  for(it in 1:2){
    for(ik in 1:10000){
      Kcal_tot_artk[1,ir,it,ik]=Kcal_tot_asrtk[1,1,ir,it,ik]+Kcal_tot_asrtk[1,2,ir,it,ik]
      Kcal_tot_artk[2,ir,it,ik]=sum(Kcal_tot_asrtk[2:40,1:2,ir,it,ik])
    }
  }
}

#calculate proportion of energy taken from bering
Bering_energy = apply(Kcal_tot_artk[,1,,],3,'sum')
Chukchi_energy = apply(Kcal_tot_artk[,2,,],3,'sum')
summary(Bering_energy/(Bering_energy+Chukchi_energy))

#calculate proportion of energy taken in summer vs fall
Summer_energy = apply(Kcal_tot_artk[,,1,],3,'sum')
Winter_energy = apply(Kcal_tot_artk[,,2,],3,'sum')
summary(Summer_energy/(Summer_energy+Winter_energy))


#data frame for plotting
Plot_df = data.frame(matrix(0,40000,5))
colnames(Plot_df)=c("Region","Season","Replicate","Series","tCal")
counter = 0
counter2 = 1
Reps = c(1:10000)
for(ir in 1:2){
  for(it in 1:2){
    cur_pl = counter+Reps
    if(ir==1)Plot_df[cur_pl,"Region"]="Bering"
    else Plot_df[cur_pl,"Region"]="Chukchi-Beaufort"
    if(it==1)Plot_df[cur_pl,"Season"]="Summer"
    else Plot_df[cur_pl,"Season"]="Winter"
    Plot_df[cur_pl,"Replicate"] = Reps
    Plot_df[cur_pl,"tCal"] = apply(Kcal_tot_artk[,ir,it,],2,'sum')
    Plot_df[cur_pl,"Series"]=counter
    counter=counter+10000
    counter2 = counter+1
  }
}
Plot_df[,"tCal"]=Plot_df[,"tCal"]/1000000000

tCal_plot = ggplot(Plot_df)+geom_density(aes(x=tCal),alpha=0.6,fill="gray")+
        facet_grid(Region~Season,scales="free")+
        scale_x_continuous(limits = c(0,250))+
        scale_y_continuous(limits = c(0,0.04))+
        xlab('Population energy requirements (tcal)')+ylab('Probability density')
tCal_plot

png("tCal_hists.png",width = 6, height = 6, units='in', res = 600)
tCal_plot
dev.off()

tCal_plot2 = ggplot(Plot_df,aes(x=Season,y=tCal,fill=Region))+geom_boxplot()+
  ylab('Population energy requirements (tcal)')
tCal_plot2


png("tCal_hists2.png",width = 6, height = 6, units='in', res = 600)
  tCal_plot2
dev.off()


#### Bring in composition data
load("./Data/prey_composition_g_kcal.RData")
#indices are MCMC rep, prey species, Age_class,Region,Season
Biomass_g_kpart = B #Biomass in grams
for(k in 1:10){
  for(ia in 1:2){
    for(ir in 1:2){
      for(it in 1:2){
        Biomass_g_kpart[,k,ia,ir,it]=B[,k,ia,ir,it]*Kcal_tot_artk[ia,ir,it,]
      }
    }
  }
}

#collapse to prey species (and version with region)
Biomass_g_kp = matrix(0,10000,10)
Biomass_g_kpr = array(0,dim=c(10000,10,2))
for(ik in 1:10){
  for(ia in 1:2){
    for(ir in 1:2){
      for(it in 1:2){
        Biomass_g_kp[,ik]=Biomass_g_kp[,ik]+Biomass_g_kpart[,ik,ia,ir,it]
        Biomass_g_kpr[,ik,ir]=Biomass_g_kpr[,ik,ir]+Biomass_g_kpart[,ik,ia,ir,it]
      }
    }
  }
}
save(Biomass_g_kpr,file="Biomass_for_Johanna.RData")
Biomass_g_mkrt = apply(Biomass_g_kpart,c(1,2,4,5),'sum')
Biomass_mt_krt = apply(Biomass_g_mkrt,c(2,3,4),'mean')/1000000
Biomass_mt_kp = Biomass_g_kp / 1000000  #metric tons
Biomass_mt_kpr = Biomass_g_kpr / 1000000  #metric tons

#summary stats
colMeans(Biomass_mt_kp)  #point estimates
sqrt(apply(Biomass_mt_kp,2,'var'))/colMeans(Biomass_mt_kp)  #CV

#plot
Plot_df = data.frame(matrix(0,10000*9,2))
colnames(Plot_df)=c("Species","Biomass")
Species = c("Pacific herring","Rainbow smelt","Pacific capelin","Sand lance","Arctic cod", 
            "Saffron cod","Walleye pollock","Pacific cod","Yellowfin sole","Other fish")
counter=0

Biomass_mt_kr = apply(Biomass_g_kpr,c(2,3),'mean')/100000000
Post_mean_mt_estimates = data.frame("species"=Species,bering_mt=Biomass_mt_kr[,1],chukchi_mt=Biomass_mt_kr[,2])
write.csv(Post_mean_mt_estimates,file="Post_mean_mt_estimates.csv")
Post_mean_mt_estimates_byseason = data.frame("species"=Species,
        bering_summer_mt=Biomass_mt_krt[,1,1],bering_winter_mt=Biomass_mt_krt[,1,2],
        chukchi_summer_mt=Biomass_mt_krt[,2,1],chukchi_winter_mt=Biomass_mt_krt[,2,2])
write.csv(Post_mean_mt_estimates_byseason,file="Post_mean_mt_estimates_by_season.csv")
for(isp in 1:10){
  Cur_reps = Reps+counter
  Plot_df[Cur_reps,"Biomass"]=Biomass_mt_kp[,isp]
  Plot_df[Cur_reps,"Species"]=Species[isp]
  counter=counter+10000
}

Plot_df$Species = factor(Plot_df$Species,levels=c("Pacific capelin","Saffron cod","Arctic cod","Pacific herring","Sand lance", 
                                                            "Rainbow smelt","Walleye pollock","Pacific cod","Yellowfin sole","Other fish")
)

#remove outliers for plotting
which_rem = which(Plot_df$Biomass>10000 & Plot_df$Species=="Pacific cod")
Plot_df[which_rem,"Biomass"]=NA

biomass_plot = ggplot(Plot_df)+geom_density(aes(x=Biomass/1000))+facet_wrap(vars(Species),scales="free")+
  xlab('Biomass (1000 mt)') + ylab('Probability density')

png("fish_biomass.png",width = 8, height = 8, units='in', res = 600)
  biomass_plot
dev.off()

#alternative plot with confidence intervals, different species on same plot
Plot_df$Biomass = Plot_df$Biomass/1000
Which_hi = c("Arctic cod","Pacific capelin","Pacific herring","Saffron cod")
Plot_df$Plot = "Low"
Plot_df$Plot[which(Plot_df$Species %in% Which_hi)]="High"

biomass_plot2 = ggplot(Plot_df)+
  geom_boxplot(aes(x=Species,y=Biomass,color=Species))+ylab('Biomass (1000 mt)')+
   theme(axis.text.x=element_blank(),text = element_text(size = 14))+
   scale_color_discrete(labels = c(
     "Pacific capelin","Saffron cod","Arctic cod","Pacific herring","Sand lance", 
     "Rainbow smelt","Walleye pollock","Pacific cod","Yellowfin sole",
     expression(italic("Other fish"))
   ))

png("fish_biomass2.png",width = 8, height = 8, units='in', res = 600)
  biomass_plot2
dev.off()

#calculate total biomass and CI
Tot_biomass_sample = rowSums(Biomass_mt_kp)
mean(Tot_biomass_sample)
quantile(Tot_biomass_sample,c(0.025,0.975))
