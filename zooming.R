setwd("file-path")  #e.g. setwd("C:\\Users\\Documents"), assuming the image file is saved in the location "C:\\Users\\Documents"

library(matlab)
library(jpeg)
library(metaSEM)

#Here we are generating the registered image for the `girl' image for r1=2 and r2=6.
#In this case, we generate the MSD values for (r1=2, r2=6) as in Table 2,
#and generate the first image in the second row of Figure 4.

#For any other image and and any other r1 and r2 value, we can generate the regitered images accordingly.
#For that we need to change the `im, im_zoom, r1 and r2' accordingly.

#Here the MSD values will not exactly be equal to the tabulated values. As we are running one instance here.
#Whereas, in the paper, we ran this for 10 times and took the average of the MSD values.
#Though the MSD values will be pretty close to the tabulated value.

#Run-time for one inastance is about 5-8 minutes for the `girl' image.


im<-readJPEG("girl.jpg")    ##reference input image, `girl image in this case. Change the image file name accordinly for any other image
im_zoom<-readJPEG("girl_z.jpg")  ##zoomed input image


## Registration algorithm for zoomed in real image##
w1<-20

img<- padarray(im_zoom,c(w1,w1),"symmetric","both")  # padded zoomed image (by the def. of the paper)
img_zoom<-padarray(im,c(w1,w1),"symmetric","both")    # padded reference image (by the def. of the paper) 

#set.seed(8)
m=proc.time()
stop<-FALSE



#Change the parameter values of r1 and r2 below. 
r1<-2
r2<-6

img<-img + matrix(rnorm(nrow(img)*ncol(img),0,0.01),nrow = nrow(img),ncol = ncol(img))
img_zoom<-img_zoom + matrix(rnorm(nrow(img)*ncol(img),0,0.01),nrow = nrow(img),ncol = ncol(img))

img_L2<-matrix(rep(NA,nrow(img)*ncol(img)),nrow = nrow(img),ncol = ncol(img))
img_L1<-matrix(rep(NA,nrow(img)*ncol(img)),nrow = nrow(img),ncol = ncol(img))
img_cor<-matrix(rep(NA,nrow(img)*ncol(img)),nrow = nrow(img),ncol = ncol(img))


for(i in (w1+1):(nrow(img)-w1))
{
  for(j in (w1+1):(ncol(img)-w1))
  {
    
    N1<-c()
    for(k1 in (i-r1):(i+r1))
    {
      for(l1 in (j-r1):(j+r1))
      {
        if((k1-i)^2 + (l1-j)^2 <= (r1)^2)
        {
          N1<-c(N1,img[k1,l1])
        }
      }
    }
    msd<-c()
    mad<-c()
    r<-c()
    xz<-c()
    yz<-c()
    cov1<-0
    sd1<-0
    sd2<-0
    for(m1 in (i-r2):(i+r2))
    {
      for(n1 in (j-r2):(j+r2))
      {
        if((m1-i)^2+(n1-j)^2 <= (r2)^2)
        {
          N2<-c()
          N<-0
          for(s1 in (m1-r1):(m1+r1))
          {
            for(t1 in (n1-r1):(n1+r1))
            {
              if((s1-m1)^2+(t1-n1)^2 <= (r1)^2)
              {
                N2<-c(N2,img_zoom[s1,t1])
                N<-N+1
              }
            }
          }
          cov1<-cov(N1,N2)
          sd1<-sqrt(var(N1))
          sd2<-sqrt(var(N2))
          r<-c(r,(cov1/(sd1*sd2)))
          msd<-c(msd,sum((N1-N2)^2)/N)
          mad<-c(mad,sum(abs(N1-N2))/N)
          xz<-c(xz,round(s1))
          yz<-c(yz,round(t1))
          if((sd1==0)|(sd2==0))
          {
            stop=TRUE
            print(sd1)
            print(sd2)
            print(N1)
            print(N2)
            break
          }
        }
        if(stop){break}
      }
      if(stop){break}
    }
    if(stop){break}
    r[is.infinite(r)]<-0
    index_L2<-which.min(msd)
    index_L1<-which.min(mad)
    index_cor<-which.max(r)
    img_L2[xz[index_L2]-r1,yz[index_L2]-r1]<- img[i,j]
    img_L1[xz[index_L1]-r1,yz[index_L1]-r1]<- img[i,j]
    img_cor[xz[index_cor]-r1,yz[index_cor]-r1]<- img[i,j]
  }
  if(stop){break}
}

n=proc.time()

n-m
##Removing Na'S

h<-2
for(i in (w1+r1):(nrow(img)-w1-r1+1))
{
  for(j in (w1+r1):(ncol(img)-w1-r1+1))
  {
    if(is.na(img_L1[i,j]))
    {
      s<-c()
      for( k in (i-h):(i+h))
      {
        for(l in (j-h):(j+h))
        {
          s<- c(s,img_L1[k,l]) 
        }
      }
      img_L1[i,j]<-mean(s,na.rm = TRUE)
    }
    if(is.na(img_L2[i,j]))
    {
      s1<-c()
      for( k in (i-h):(i+h))
      {
        for(l in (j-h):(j+h))
        {
          s1<- c(s1,img_L2[k,l]) 
        }
      }
      img_L2[i,j]<-mean(s1,na.rm = TRUE)
    }
    if(is.na(img_cor[i,j]))
    {
      s2<-c()
      for( k in (i-h):(i+h))
      {
        for(l in (j-h):(j+h))
        {
          s2<- c(s2,img_cor[k,l]) 
        }
      }
      img_cor[i,j]<-mean(s2,na.rm = TRUE)
    }
  }
}

#MSE 
t1<-0
t2<-0
t3<-0
n1<-0
n2<-0
n3<-0
for(i in (w1+r1):(nrow(img)-w1-r1+1))
{
  for(j in (w1+r1):(nrow(img)-w1-r1+1))
  {
    if(is.na(img_L1[i,j])==FALSE)
    {
      t1<-t1+(img_zoom[i,j]-img_L1[i,j])^2
      n1<-n1+1
    }
    
    if(is.na(img_L2[i,j])==FALSE)
    {
      t2<-t2+(img_zoom[i,j]-img_L2[i,j])^2
      n2<-n2+1
    }
    
    if(is.na(img_cor[i,j])==FALSE)
    {
      t3<-t3+(img_zoom[i,j]-img_cor[i,j])^2
      n3<-n3+1
    }
    
    
  }
}

MSD_L1<-sqrt(t1/n1)
MSD_L2<-sqrt(t2/n2)
MSD_cor<-sqrt(t3/n3)


L1_img<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)  #Registered image under L1-norm
L1_img[(r1+1):((nrow(L1_img)-r1)),(r1+1):((ncol(L1_img)-r1))]<-img_L1[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

L2_img<--matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1) #Registered image under L2-norm
L2_img[(r1+1):((nrow(L2_img)-r1)),(r1+1):((ncol(L2_img)-r1))]<-img_L2[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

cor_img<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1) #Registered image under CC-method
cor_img[(r1+1):((nrow(cor_img)-r1)),(r1+1):((ncol(cor_img)-r1))]<-img_cor[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]




##Anomaly detection

anomaly_L2<-matrix(data=NA,nrow = nrow(img),ncol = ncol(img))
anomaly_L1<-matrix(data=NA,nrow = nrow(img),ncol = ncol(img))
anomaly_cor<-matrix(data=NA,nrow = nrow(img),ncol = ncol(img))

for(i in (w1+r1+1):(nrow(img)-w1-r1))
{
  for(j in (w1+r1+1):(ncol(img)-w1-r1))
  {
    if((abs(img_zoom[i,j]-img_L1[i,j])< (max(max(img))-min(min(img)))*0.15))
      anomaly_L1[i,j]<-0
    
    else
      anomaly_L1[i,j]<-1
    
    if((abs(img_zoom[i,j]-img_L2[i,j])< (max(max(img))-min(min(img)))*0.15))
      anomaly_L2[i,j]<-0
    
    else
      anomaly_L2[i,j]<-1
    
    if((abs(img_zoom[i,j]-img_cor[i,j])< (max(max(img))-min(min(img)))*0.15))
      anomaly_cor[i,j]<-0
    
    else
      anomaly_cor[i,j]<-1
    
  }
}  


L1_anmy<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)  #Anomaly image under L1-norm
L1_anmy[(r1+1):((nrow(L1_anmy)-r1)),(r1+1):((ncol(L1_anmy)-r1))]<-anomaly_L1[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

L2_anmy<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)   #Anomaly image under L2-norm
L2_anmy[(r1+1):((nrow(L2_anmy)-r1)),(r1+1):((ncol(L2_anmy)-r1))]<-anomaly_L2[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

cor_anmy<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)   #Anomaly image under CC-method
cor_anmy[(r1+1):((nrow(cor_anmy)-r1)),(r1+1):((ncol(cor_anmy)-r1))]<-anomaly_cor[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

windows(10,10)
image(rot90(L1_img,3),col = grey(seq(0,1,length=256))) ##Registered image for L1-norm.

MSD_L1
MSD_L2
MSD_cor
