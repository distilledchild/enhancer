
library(plyr)
library(ggplot2)
library(broom)
library(chron)
library(readr)
library(purrr)
library(stringr)
library(dplyr)
library(stringr)



#function for quantile normalizing data
normf<- function(x){
  idx<- !is.na(x)
  o<- order(x[idx])
  p<- (1:sum(idx)-0.5)/sum(idx)
  z<- qnorm(p)
  x[idx]<- z[order(o)]
  x
}





raw=read.csv("data/raw_data/raw_data.csv",header=T,stringsAsFactors = F,na.strings=c("","NA"))





#convert rfids to character
raw$rfid <- as.character(raw$rfid)


#load traits and covariates
load("data/raw_data/traits.RData")



#raw<-raw %>%
#  mutate_if(is.character, str_trim)

if(length(grep("project_name",colnames(raw)))==1){
  meta_analysis <- TRUE
}



centers<-unique(raw$project_name)
sexes<-unique(raw$sex)

qnormed_center_data<-list()

for(l in 1:length(centers)){

  data<-raw[which(raw$project_name==centers[l]),]


  qnormed_sex_data_all <- list()

  for(j in 1:length(sexes)){

    qnormed_data_each_trait<- list()
    for(i in 1:length(traits)){

      temp<-data[which(data$sex==sexes[j]),c("rfid",traits[i])]

      #may be try a regex for white spaces
      #this part removes dirty data : white spaces NANs etc
      temp<-temp[which(!temp[,traits[i]]==""),]
      temp<-temp[which(!(is.infinite(temp[,traits[i]]) | is.na(temp[,traits[i]]) |is.nan(temp[,traits[i]]))),c("rfid",traits[i])]

      temp[,paste0("qnormed_",traits[i])] <-normf(temp[,traits[i]])


      qnormed_data_each_trait[[i]]<-temp[,c("rfid",paste0("qnormed_",traits[i]))]

    }

    qnormed_sex_data_all[[j]] <-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "rfid", all = TRUE),
                                       qnormed_data_each_trait)

  }

  qnormed_sex_data_all<-do.call("rbind",qnormed_sex_data_all)
  qnormed_center_data[[l]]<-merge(qnormed_sex_data_all,data[,c("rfid",covs)],by="rfid",all=T)
  names(qnormed_center_data)[l]<-centers[l]
}



rm(qnormed_data_each_trait,qnormed_sex_data_all)

#check for PVE by age




qnormed_traits<-paste0("qnormed_",traits)


age_covs<-covs[grepl("_age",covs)]
summary_age_all <- list()

for(l in 1:length(centers)){


  summary_age_each_covariate<-list()
  data<-qnormed_center_data[[centers[l]]]

  for(i in 1:length(age_covs)){


    exp_var<-gsub('[[:digit:]]+', '', strsplit(age_covs[i],split="_age")[[1]][1])
    traits_for_age_var<-grep(paste0("_",exp_var),qnormed_traits,value = T)
    age_var<-age_covs[i]

    temp<-data[,c("rfid",c(age_covs[i],traits_for_age_var))]
    summary_age_each_trait <- list()
    for(j in 1:length(traits_for_age_var)){



      temp[,"age"]<-temp[,age_var]

      #remove missing age values
      temp<-temp[which(!temp[,"age"]==""),]
      temp<-temp[which(!(is.infinite(temp[,"age"]) | is.na(temp[,"age"]) |is.nan(temp[,"age"]))),]



      model<-glance(lm(temp[,traits_for_age_var[j]] ~ age, data=temp))
      df <- data.frame(model)
      df$trait=traits_for_age_var[j]
      df$age_cov=age_var
      summary_age_each_trait[[j]]<-df

    }
    summary_age_each_covariate[[i]]<-do.call("rbind",summary_age_each_trait)

  }

  all_age_summary<-do.call("rbind",summary_age_each_covariate)
  all_age_summary$r.squared<-all_age_summary$r.squared *100

  summary_age_all[[l]]<-all_age_summary[which(all_age_summary$r.squared > 2 & all_age_summary$p.value < 0.001),]
}



rm(summary_age_each_covariate,summary_age_each_trait,df,data,model)
for(k in 1:length(centers)){
  summary_age_all[[k]]$center<-centers[k]
}

summary_lm_age_to_regress<-summary_age_all
summary_lm_age_to_regress<-do.call("rbind",summary_lm_age_to_regress)


save(summary_lm_age_to_regress,file="pheno_processing_summary/covs/age_regressed.RData")


names(summary_age_all)<-centers


for(l in 1:length(centers)){


  if(nrow(summary_age_all[[l]])>0){

    age_regressed_out_temp<-summary_age_all[[centers[l]]]

    data<-qnormed_center_data[[centers[l]]]

    temp_list<-list()
    for(i in 1:nrow(age_regressed_out_temp)){
      temp<-data
      age_var<-age_regressed_out_temp$age_cov[i]
      traits_age<-age_regressed_out_temp$trait[i]
      temp[,"age"]<-temp[,age_var]

      temp_list[[i]]<- temp[,c("rfid","age",traits_age)]

      temp1<-temp_list[[i]]

      temp1<-temp1[which(!(is.infinite(temp1[,traits_age]) | is.na(temp1[,traits_age]) |is.na(temp1[,"age"])| is.nan(temp1[,traits_age]))),c("rfid","age",traits_age)]
      temp_list[[i]]<-temp1


      temp_list[[i]][,paste0("resid_age_",traits_age)]<- resid(lm(get(traits_age) ~ age, temp_list[[i]]))

      #qnorm the residuals
      temp_list[[i]][,paste0("qnormed_resid_age_",traits_age)] <-normf(temp_list[[i]][,paste0("resid_age_",traits_age)])
      data[,traits_age]<-NULL
      data<-merge(data,temp_list[[i]][,c("rfid",paste0("qnormed_resid_age_",traits_age))],by="rfid",all=T)




    }
    qnormed_center_data[[centers[l]]]<-data
  }

  colnames(qnormed_center_data[[centers[l]]])



}






#clean up column names


for(l in 1:length(centers)){

  temp<-qnormed_center_data[[centers[l]]]
  colnames(temp)<-gsub("qnormed_resid_age_qnormed_","",colnames(temp))
  colnames(temp)<-gsub("qnormed_","",colnames(temp))
  qnormed_center_data[[centers[l]]]<-temp

}



#find discrete covariates


#check and regress out the effect of cohort
#we use a separate code chunk for cohort because it is common across all experiments

discrete_covs<-names(grep("character",sapply(covs, function(x) class(raw[[x]])),value=T) )

#remove sex because we've quantile normalized the sexes separately. This is similar to accounting
#for sex as a covariate


discrete_covs<-grep("sex",discrete_covs,invert = T,value = T)
discrete_covs<-grep("cohort",discrete_covs,invert = T,value = T)

covs_common <- c("cohort")



if(length(covs_common)>0){

  center_cov_traits<-list()
  for(l in 1:length(centers)){


    data <- qnormed_center_data[[l]]

    #remove white spaces etc from cohort column values
    center_cov_data<-list()
    for(i in 1:length(covs_common)){

      temp_df<-data

      temp_df<-temp_df[,c("rfid",covs_common[i],traits)]


      temp_df<-temp_df[complete.cases(temp_df[ , c("rfid",covs_common[i])]), ]

      df_name<-paste0("forBinary_",covs_common[i])
      assign("temp_covs",data.frame(temp_df[,covs_common[i]]))
      colnames(temp_covs)<-covs_common[i]
      temp_cat = data.frame (model.matrix (~ temp_covs[,covs_common[i]]-1, data=temp_df))
      names(temp_cat)<-levels(temp_covs[,covs_common[i]])
      assign(paste0(covs_common[i],"_cat"),temp_cat)
      temp_df <- cbind (temp_df, temp_cat)




      master <- list()
      master <- vector("list",length(traits))
      names(master)<- traits



      master_summary<-setNames(replicate(length(traits),data.frame()),traits)






      for ( k in 1:length(traits)){

        models<-list()

        cnames<- unique(temp_df[,covs_common[i]])

        models <- lapply (cnames, function(x){
          anova (lm (get(traits[k]) ~ 1,data=temp_df), lm (substitute (get(traits[k]) ~ r, list (r= as.name(x))), data=temp_df) )
        })
        master[[k]] <- models


        for (n in 1:length(cnames)){

          p_val<-master[[i]][[n]]$Pr[2]
          pve <- ((master[[i]][[n]]$RSS[1] - master[[i]][[n]]$RSS[2])/master[[i]][[n]]$RSS[1]) * 100


          df <- data.frame(p_val,pve)

          master_summary[[k]] <- rbind (master_summary[[k]],df)
        }
        master_summary[[k]]$cnames <- cnames
        master_summary[[k]]$trait<- rep(traits[k],length(cnames))



      }
      #master_summary[[k]]<-master_summary[[k]][which(master_summary[[k]]$pve > 2 & master_summary[[k]]$p_val < 0.001),]

      # master_summary<-do.call("rbind",master_summary)
      center_cov_data[[i]]<-master_summary



    }
    center_cov_traits[[l]]<-center_cov_data
  }
  need_to_regress<-list()

  for(l in 1:length(centers)){
    data<-center_cov_traits[[l]]

    test<-do.call(Map, c(f = rbind,data))
    final_summary_covs<-do.call("rbind",test)
    need_to_regress[[l]]<-final_summary_covs[which(final_summary_covs$pve >2 & final_summary_covs$p_val < 0.001),]
    names(need_to_regress)[l]<-centers[l]
  }

  need_to_regress<-do.call("rbind",need_to_regress)
  write.csv(need_to_regress,"pheno_processing_summary/covs/common_covariates_regressed.csv",row.names = F,quote=T)





}






#check and regress out the effect of other covariates. These are specific to each experiment

if(length(discrete_covs)>0){


  center_cov_traits<-list()


  for(l in 1:length(centers)){

    data <- qnormed_center_data[[l]]

    #remove white spaces etc from cohort column values
    center_cov_data<-list()

    for(i in 1:length(discrete_covs)){

      temp_df<-qnormed_center_data[[l]]

      temp_df<-temp_df[,c("rfid",discrete_covs[i],traits)]


      temp_df<-temp_df[complete.cases(temp_df[ , c("rfid",discrete_covs[i])]), ]

      df_name<-paste0("forBinary_",discrete_covs[i])
      assign("temp_covs",data.frame(temp_df[,discrete_covs[i]]))
      colnames(temp_covs)<-discrete_covs[i]


      temp_cat = data.frame (model.matrix (~ temp_covs[,discrete_covs[i]]-1, data=temp_df))
      names(temp_cat)<-levels(temp_covs[,discrete_covs[i]])
      assign(paste0(discrete_covs[i],"_cat"),temp_cat)
      temp_df <- cbind (temp_df, temp_cat)

      exp_var<-gsub('[[:digit:]]+', '', strsplit(discrete_covs[i],split="_")[[1]][1])
      traits_for_covs<-grep(paste0(exp_var,"_"),traits,value = T)

      master <- list()
      master <- vector("list",length(traits_for_covs))
      names(master)<- traits_for_covs



      master_summary<-setNames(replicate(length(traits_for_covs),data.frame()),traits_for_covs)






      for ( k in 1:length(traits_for_covs)){

        models<-list()

        cnames<- unique(temp_df[,discrete_covs[i]])

        models <- lapply (cnames, function(x){
          anova (lm (get(traits_for_covs[k]) ~ 1,data=temp_df), lm (substitute (get(traits_for_covs[k]) ~ r, list (r= as.name(x))), data=temp_df) )
        })
        master[[k]] <- models


        for (n in 1:length(cnames)){

          p_val<-master[[i]][[n]]$Pr[2]
          pve <- ((master[[i]][[n]]$RSS[1] - master[[i]][[n]]$RSS[2])/master[[i]][[n]]$RSS[1]) * 100


          df <- data.frame(p_val,pve)

          master_summary[[k]] <- rbind (master_summary[[k]],df)
        }
        if(nrow(master_summary[[k]])>0){
          master_summary[[k]]$cnames <- cnames
          master_summary[[k]]$trait<- rep(traits_for_covs[k],length(cnames))
          master_summary[[k]]$cov<- rep(discrete_covs[i],length(cnames))
        }




      }

      center_cov_data[[i]]<-master_summary



    }
    center_cov_traits[[l]]<-center_cov_data
  }

}





need_to_regress<-list()

rm(center_cov_data)

for(l in 1:length(centers)){
  data<-center_cov_traits[[l]]

  test<-do.call(Map, c(f = rbind,data))
  final_summary_covs<-do.call("rbind",test)
  need_to_regress[[l]]<-final_summary_covs[which(final_summary_covs$pve >2 & final_summary_covs$p_val < 0.001),]
  names(need_to_regress)[l]<-centers[l]
}




if(length(need_to_regress)>0){


  for(l in 1:length(centers)){

    center_specific<-need_to_regress[[centers[l]]]
    center_specific<-unique(center_specific)
    counts<-plyr::count(center_specific,vars="trait")
    center_specific<-center_specific[order(center_specific$trait),]


    for (i in seq_along(counts$trait)){
      filename<-paste0("lmformula.R")
      if (file.exists(filename))
        #Delete file if it exist
        file.remove(filename)


      idx<-match(counts$trait[i],center_specific$trait)


      f<- counts$freq[i]
      last<- (idx) + (f-1)
      center_specific[idx:last,]
      cvs<-center_specific$cnames[idx:last]
      cov<-center_specific$cov[idx:last]

      temp_df<-qnormed_center_data[[centers[l]]]



      temp_df<-temp_df[,c("rfid",cov,counts$trait[i])]


      temp_df<-temp_df[complete.cases(temp_df[ , c("rfid",cov)]), ]

      df_name<-paste0("forBinary_",cov)
      assign("temp_covs",data.frame(temp_df[,cov]))
      colnames(temp_covs)<-cov


      temp_cat = data.frame (model.matrix (~ temp_covs[,cov]-1, data=temp_df))
      names(temp_cat)<-levels(temp_covs[,cov])
      assign(paste0(cov,"_cat"),temp_cat)
      temp_df <- cbind (temp_df, temp_cat)
      temp_df<-temp_df[complete.cases(temp_df), ]

      idx<-match(counts$trait[i],center_specific$trait)
      part1 <- paste0("temp_df$resid_",counts$trait[i],"= resid( lm ( ",counts$trait[i]," ~ ",center_specific$cnames[idx]," ")


      f<- counts$freq[i]
      if (f > 1){
        after<-idx+1
        last<- (idx) + (f-1)
        covs<-center_specific$cnames[after: last]
        plus_after <- "+"
        with_plus<-paste(covs,collapse="+")
        formula <- paste0 (part1,plus_after,with_plus,",data=temp_df)")
        write(formula,file=filename,append=TRUE)


      }
      else{
        formula <- paste0(part1,",data=temp_df))")
        write(formula,file=filename,append=TRUE)
      }
      source(file =filename)
      v<-paste0("norm_",counts$trait[i])

      resid<-paste0("resid_",counts$trait[i])
      temp_df[,v] <-normf(temp_df[,resid])
      temp_df<-temp_df[,c("rfid",v)]
      qnormed_center_data[[centers[l]]][,counts$trait[i]]<-NULL
      qnormed_center_data[[centers[l]]]<-merge(temp_df,qnormed_center_data[[centers[l]]],by="rfid",all=T)
    }




  }

}





colnames(qnormed_center_data[[1]])<-gsub("norm_","",colnames(qnormed_center_data[[1]]))
colnames(qnormed_center_data[[2]])<-gsub("norm_","",colnames(qnormed_center_data[[2]]))

colnames(qnormed_center_data[[1]])<-gsub("qnorm_","",colnames(qnormed_center_data[[1]]))
colnames(qnormed_center_data[[2]])<-gsub("qnorm_","",colnames(qnormed_center_data[[2]]))



qnormed_center_data<-do.call("rbind",qnormed_center_data)



gwas<-qnormed_center_data


save(gwas,file="data/residuals/residuals.RData")
