 
remove(list = ls())


setwd("C:\\Users\\cssc\\OneDrive\\Documents\\Registration for Zooming") 
  library(matlab)
  library(jpeg)
  ##Construction of original image
  r<- 64
  c<- 64
  sd<-0.01
  x1<-seq(0,1,length=r)
  y<-seq(0,1,length=c)
  el_b<-rep(0,r^2)
  m_b<-matrix(el_b,nrow = r,byrow = TRUE)  #elements of black region
  el_w<-rep(1,(round(3*r/4)-round(r/4))^2) #elements of white region
  m_w<-matrix(el_w,nrow = round(3*r/4)-round(r/4),byrow = TRUE)
  m_b[(round(r/4)+1):round(3*r/4),(round(r/4)+1):round(3*r/4)]<-m_w
  
 
  
  dim(m_w)
  
  windows(10,10)
  image(m_b,col = grey(seq(0,1,length=256)))
  
  ##Construction of zoomed image
  
  # zoom_low<-m_b[7:58,7:58]
  # 
  # jpeg(filename = "E:\\zoom.jpg",height=64,width=64)
  # par(bty="n",xaxt="n",yaxt="n",mar=c(0,0,0,0))
  # image(zoom_low,col = grey(seq(0,1,length=256)))
  # dev.off()
  # 
  # zoom<-readJPEG("zoom.jpg")
  
  x=matrix(0,64,64)
  x[11:54,11:54]=1
  m_b=x+matrix(rnorm(r*c,0,sd),64,64)
  
  y=matrix(0,64,64)
  y[6:59,6:59]=1
  m_b_zoom=y+matrix(rnorm(r*c,0,sd),64,64)
  
  par(mfrow=c(1,2))
  windows(10,10)
  image(m_b_zoom)
  
  ##Adding noise to both images
  m_b_n<-matrix(NA,nrow = nrow(m_b),ncol = ncol(m_b))
  m_b_n<-m_b+ matrix(rnorm(nrow(m_b)^2,0,sd),nrow = nrow(m_b),ncol = ncol(m_b))
  
  
  ##Work on girl image
  
  girl<-readJPEG("girl.jpg")
  dim(girl)
  
  girl_z<-girl[7:122,7:122]
  
  jpeg(filename = "stop__zoom_noise.jpg",height=128,width=128)
   par(bty="n",xaxt="n",yaxt="n",mar=c(0,0,0,0))
   image(rot90(stop_z_n,3),col = grey(seq(0,1,length=256)))
   dev.off()
   
   girl_zm<-readJPEG("girl_z.jpg")
  
   ##Adding noise
   pz<-stop_z[,,1]
   
  stop_n<-stop[,,1]+matrix(rnorm(nrow(stop[,,1])*ncol(stop[,,1]),0,0.05),nrow(stop[,,1]),ncol(stop[,,1]))
   stop_z_n<-pz+matrix(rnorm(nrow(pz)*ncol(pz),0,0.05),nrow(pz),ncol(pz))
   dim(stop_z)
  
  windows(10,10)
  par(mfrow=c(1,2))
  image(m_b_n,col = grey(seq(0,1,length=256)))
  image(zoom,col = grey(seq(0,1,length=256)))
  
  
  ##Moon image
  moon<-  readJPEG("moon.jpg")
  moon_l<-moon[7:122,7:122]
  jpeg(filename = "moon_z.jpg",height=128,width=128)
  par(bty="n",xaxt="n",yaxt="n",mar=c(0,0,0,0))
  image(rot90(moon_l,3),col = grey(seq(0,1,length=256)))
  dev.off()
  
  
  #Mountain image
  mountain<-  readJPEG("mountain.jpg")
  mountain_l<-mountain[7:122,7:122]
  jpeg(filename = "girl_z.jpg",height=128,width=128)
  windows(10,10)
  par(bty="n",xaxt="n",yaxt="n",mar=c(0,0,0,0))
  image(rot90(girl_zm,3),col = grey(seq(0,1,length=256)))
  box(which = "figure")
  dev.off()
  
  mountain_z<-readJPEG("mountain_z.jpg")
  
  
  ##Stop image
  stop<-readJPEG("stop.jpg")
  stop_l<-stop[7:122,7:122,1]
  jpeg(filename = "stop_z.jpg",height=128,width=128)
  par(bty="n",xaxt="n",yaxt="n",mar=c(0,0,0,0))
  image(rot90(stop_l,3),col = grey(seq(0,1,length=256)))
  dev.off()
  stop_z<-readJPEG("stop_z.jpg")
  
  
  ##Work on pepper image
  
  pepper<-readJPEG("test_pepper.jpg")
  
  dim(pepper)
  
  pepper_z<-pepper[6:251,6:251]
  
  jpeg(filename = "pepper_lz.jpg",height=256,width=256)
  par(bty="n",xaxt="n",yaxt="n",mar=c(0,0,0,0))
  image(rot90(pepper_z,3),col = grey(seq(0,1,length=256)))
  dev.off()
  
  pepper_zm<-readJPEG("pepper_lz.jpg")
  dim(pepper_zm)
  
  par(mfrow=c(1,2))
  image(rot90(pepper,3),col = grey(seq(0,1,length=256)))
  image(rot90(pepper_zm,3),col = grey(seq(0,1,length=256)))
  
  
  
  
  ##Finding a initial estimate for r2
  ref_edge<-which(stepEdgeLC2K(image=m_b_n,bandwidth=3,thresh=0.2,plot=TRUE)==1,arr.ind = TRUE)
  
  zoom_edge<- which(stepEdgeLC2K(image=zoom,bandwidth=3,thresh=0.2,plot=TRUE)==1,arr.ind = TRUE)
  
 
  
  adec<-function(mat1,c1,c2)  ##function to calculate average distance of edges from center
  {
  s<-c()
  
  for(i in 1:nrow(mat1))
  {
    s<-c(s,sqrt((mat1[i,1]-c1)^2+(mat1[i,2]-c2)^2))
  }
  return(c(max(s),min(s)))

  }
  
  avg_ref<-adec(ref_edge,32,32)
  
  avg_zoom<-adec(zoom_edge,32,32)
  
  dif<- (abs(avg_zoom[2] - avg_ref[2])+abs(avg_zoom[1] - avg_ref[1]))/2
 ############################################################
  adec<-function(mat1,c1,c2)  ##function to calculate average distance of edges from center
  {
    s<-0
    n<-0
    for(i in 1:nrow(mat1))
    {
      if((mat1[i,1] < c1) && (mat1[i,2]>c2))
      {
      s<-s+sqrt((mat1[i,1]-c1)^2+(mat1[i,2]-c2)^2)
      n<- n+1
      }
    }
    return(s/n)
    
  }
  
  avg_ref<-adec(ref_edge,32,32)
  
  avg_zoom<-adec(zoom_edge,32,32)
  
  dif<- avg_zoom - avg_ref
###############################################
  ##Simulated image for zooming algorithm
  
  x=matrix(0,20,20)
  x[8:13,8:13]=1
  m_b=x+matrix(rnorm(20*20,0,0.01),20,20)
  windows(10,10)
  image(rot90(m_b,3),col=grey(seq(0,1,length=256)))
  
  y=matrix(0,20,20)
  y[6:15,6:15]=1
  m_b_zoom=y+matrix(rnorm(20*20,0,0.01),20,20)
  windows(10,10)
  image(rot90(m_b_zoom,3),col=grey(seq(0,1,length=256)))
###############################################  
    
  ##Finding edge pixel
  
  # library(matlab)
  # m_b_pad<-padarray(m_b,c(15,15),"symmetric","both")  ##Creation of padded image
  # 
  # w1=15
  # 
  # n<-matrix(rep(0,128*128),nrow = 128)
  # a<-matrix(rep(0,128*128),nrow = 128)
  # for(i in (w1+1):(nrow(m_b_pad)-w1))
  # {
  #   for(j in (w1+1):(nrow(m_b_pad)-w1))
  #   {
  #     for(k in (i-w1):(i+w1))
  #     {
  #       for(l in (j-w1):(j+w1))
  #       {
  #         if((i-k)^2 + (j-l)^2 < (w1^2))
  #         {
  #           n[i-w1,j-w1]<- n[i-w1,j-w1]+m_b_pad[k,l]
  #           a[i-w1,j-w1]<- a[i-w1,j-w1]+1
  #         }
  #       }
  #     }
  #   }
  # }
  # 
  # x<-n/a
  # 
  # library(metaSEM) ##contains function "vec2symMat"
  # 
  # 
  # 
  # b_0<-matrix(NA,(nrow(m_b_pad)-2*w1),(ncol(m_b_pad)-2*w1))
  # b_1<-matrix(NA,(nrow(m_b_pad)-2*w1),(ncol(m_b_pad)-2*w1))
  # b_2<-matrix(NA,(nrow(m_b_pad)-2*w1),(ncol(m_b_pad)-2*w1))
  # 
  # for(i in (w1+1):(nrow(m_b_pad)-w1))
  # {
  #   for(j in (w1+1):(nrow(m_b_pad)-w1))
  #   {
  #     
  #     M11<-0
  #     M12<-0
  #     M13<-0
  #     M22<-0
  #     M23<-0
  #     M33<-0
  #     v1<-0
  #     v2<-0
  #     v3<-0
  #     M<-NULL
  #     for(k in (i-w1):(i+w1))
  #     {
  #       for(l in (j-w1):(j+w1))
  #       {
  #         if((i-k)^2 + (j-l)^2 < (w1^2))
  #         {
  #           M11<-M11+1
  #           M12<-M12+(k-i)
  #           M13<-M13+(l-j)
  #           M22<-M22+(k-i)^2
  #           M23<-M23+(k-i)*(l-j)
  #           M33<-M33+(l-j)^2
  #           
  #           v1<-v1+m_b_pad[k,l]
  #           v2<-v2+m_b_pad[k,l]*(k-i)
  #           v3<-v3+m_b_pad[k,l]*(l-j)
  #         }
  #         
  #       }
  #     }
  #     M<-vec2symMat(c(M11,M12,M13,M22,M23,M33))
  #     b_0[i-w1,j-w1]<-(solve(M) %*% c(v1,v2,v3))[1]
  #     b_1[i-w1,j-w1]<-(solve(M) %*% c(v1,v2,v3))[2]
  #     b_2[i-w1,j-w1]<-(solve(M) %*% c(v1,v2,v3))[3]
  #   }
  # }
  # 
  # 
  # w2<-2
  # 
  # 
  # 
  # S<-matrix(rep(0,128*128),nrow = nrow(m_b_pad)-2*w1,ncol = ncol(m_b_pad)-2*w1)
  # S_N1<-matrix(rep(0,128*128),nrow = nrow(m_b_pad)-2*w1,ncol = ncol(m_b_pad)-2*w1)
  # S_N2<-matrix(rep(0,128*128),nrow = nrow(m_b_pad)-2*w1,ncol = ncol(m_b_pad)-2*w1)
  # pts<-matrix(rep(0,128*128),nrow = nrow(m_b_pad)-2*w1,ncol = ncol(m_b_pad)-2*w1)
  # pts_N1<-matrix(rep(0,128*128),nrow = nrow(m_b_pad)-2*w1,ncol = ncol(m_b_pad)-2*w1)
  # pts_N2<-matrix(rep(0,128*128),nrow = nrow(m_b_pad)-2*w1,ncol = ncol(m_b_pad)-2*w1)
  # 
  # 
  # 
  # 
  # for(i in (w1+1):(nrow(m_b_pad)-w1))
  # {
  #   for(j in (w1+1):(nrow(m_b_pad)-w1))
  #   {
  #     theta<-atan(b_2[i-w1,j-w1]/b_1[i-w1,j-w1])
  #     px_N1<-  as.integer(round(i+2*w2*cos(theta),digits = 0))
  #     py_N1<- as.integer(round(j+2*w2*sin(theta),digits = 0))
  #     px_N2<- as.integer(round(i-2*w2*cos(theta),digits = 0))
  #     py_N2<- as.integer(round(j-2*w2*sin(theta),digits = 0))
  #     
  #     
  #     for(k in (i-w2):(i+w2))
  #     {
  #       for (l in (j-w2):(j+w2))
  #       {
  #         if((i-k)^2+(j-l)^2< (w2)^2)
  #         {
  #           S[i-w1,j-w1]<-S[i-w1,j-w1]+m_b_pad[k,l]
  #           pts[i-w1,j-w1]<-pts[i-w1,j-w1]+1
  #         }
  #       }
  #     }
  #     
  #     
  #     for(k in (px_N1-w2):(px_N1+w2))
  #     {
  #       for (l in (py_N1-w2):(py_N1+w2))
  #       {
  #         if((px_N1-k)^2+(py_N1-l)^2< (w2)^2)
  #         {
  #           S_N1[i-w1,j-w1]<-S_N1[i-w1,j-w1]+m_b_pad[k,l]
  #           pts_N1[i-w1,j-w1]<-pts_N1[i-w1,j-w1]+1
  #         }
  #       }
  #     }
  #     
  #     for(k in (px_N2-w2):(px_N2+w2))
  #     {
  #       for (l in (py_N2-w2):(py_N2+w2))
  #       {
  #         if((px_N2-k)^2+(py_N2-l)^2< (w2)^2)
  #         {
  #           S_N2[i-w1,j-w1]<-S_N2[i-w1,j-w1]+m_b_pad[k,l]
  #           pts_N2[i-w1,j-w1]<-pts_N2[i-w1,j-w1]+1
  #         }
  #       }
  #     }
  #     
  #   }
  # }
  # 
  # 
  # avg_S<-S/pts
  # avg_N1<-S_N1/pts_N1
  # avg_N2<-S_N2/pts_N2
  # 
  # 
  # sigma<- sqrt(mean(as.vector((m_b-x)^2)))
  # thres<-(sqrt(2)*sigma*qnorm((1-(0.001/2)),0,1))/(2*w2+1)
  # 
  # delta<-matrix(data=NA,nrow=nrow(m_b),ncol=nrow(m_b))
  # 
  # jump<-NULL
  # 
  # 
  # edge<-matrix(data=NA,nrow=128,ncol=128)
  # 
  # for(i in 1:nrow(m_b))
  # {
  #   for(j in 1: nrow(m_b))
  #   {
  #     delta[i,j]<-min(abs(avg_S[i,j]-avg_N1[i,j]),abs(avg_S[i,j]-avg_N2[i,j]))
  #     if(delta[i,j]>thres)
  #       edge[i,j]<-1
  #     else
  #       edge[i,j]<-0
  #   }
  # }
  # 
  # windows(10,10)
  # image(edge,col = grey(seq(0,1,length=256)))
  
  
  
  
  #r1=3
  #r2=2*r1
  
  # for(i in (2*r1+r2+1):(nrow(m_b)-2*r1-r2))
  # {
  #   for(j in (2*r1+r2+1):(nrow(m_b)-2*r1-r2))
  #   {
  #       if(edge[i,j]==1)
  #       {
  #       for(k in (i-r1):(i+r1))
  #       {
  #         for(l in (j-r1):(j+r1))
  #         {
  #           if((k-i)^2 + (l-j)^2 < (r1)^2)
  #           {
  #             for(m in (k-r2):(k+r2))
  #             {
  #               for(n in (l-r2):(l+r2))
  #               {
  #                 if((m-k)^2+(n-l)^2<(r2)^2)
  #                 {
  #                   S<-0
  #                   N<-0
  #                   MSD<-c()
  #                   x1<-c()
  #                   y1<-c()
  #                   x_z<-c()
  #                   y_z<-c()
  #                   index<-0
  #                   for(s in (m-r1):(m+r1))
  #                   {
  #                     for(t in (n-r1):(n+r1))
  #                     {
  #                      if((s-m)^2+(t-n)^2 <(r1)^2)
  #                      {
  #                         S<- S+ (m_b_zoom[s,t]-m_b[s-m+i,t-n+j])^2
  #                         N<- N+1
  #                         MSD<-c(MSD,S/N)
  #                         x1<-c(x1,i)
  #                         y1<-c(y1,j)
  #                         x_z<-c(x_z,s)
  #                         y_z<-c(y_z,t)
  #                      }
  #                    }
  #                   }
  #                   index<-which.min(MSD)
  #                   m_b_zoom_reg[x1[index],y1[index]]<- m_b_zoom[x_z[index],y_z[index]]
  #                }
  #              }
  #            }
  #          }
  #        }
  #       }
  #       }
  #     else
  #     {
  #       m_b_zoom_reg[i,j]<-m_b[i,j]
  #     }
  #   }
  # }
  
  # windows(10,10)
  # par(mfrow=c(1,3))
  # image(m_b,col=gray(seq(0,1,length=256)))
  # image(m_b_zoom,col=gray(seq(0,1,length=256)))
  # image(m_b_zoom_reg,col = grey(seq(0,1,length=256)))
  
  
  ## Registration algorithm for zoomed in image ##
  
  library(metaSEM)
  w1<-25  ##pad-width
  m_b_pad<-  padarray(m_b,c(w1,w1),"symmetric","both")         #lake_trz_pad +rnorm(158*158,0,0.01)
  zoom_pad<-padarray(m_b_zoom,c(w1,w1),"symmetric","both")         # lake_ref_pad +rnorm(158*158,0,0.01)
  
  dim(m_b_zoom)
  
  girl<-readJPEG("E:\\girl.jpg")
  dim(girl)
  
  grl_z<-girl[9:120,9:120]
  
 
  
  jpeg(filename = "E:\\grl_z.jpg",height=128,width=128)
  par(bty="n",xaxt="n",yaxt="n",mar=c(0,0,0,0))
  image(rot90(grl_z,3),col = grey(seq(0,1,length=256)))
  dev.off()
  
  girl_z<-readJPEG("E:\\grl_z.jpg")
  
  windows(10,10)
  image(rot90(girl_z,3),col = grey(seq(0,1,length=256)))
  dim(girl)
  
  girl_pad <-  padarray(girl_n,c(w1,w1),"symmetric","both") #+rnorm((nrow(girl)+2*w1)*(nrow(girl)+2*w1),0,0.01)
  girl_z_pad<-padarray(girl_zm_n,c(w1,w1),"symmetric","both")#+rnorm((nrow(girl)+2*w1)*(nrow(girl)+2*w1),0,0.01)
  
  mat1<-girl_z_pad
  mat2<-girl_pad
  
  m_b_zoom_reg_L2<-matrix(data = NA,nrow = nrow(mat1),ncol = ncol(mat1))
  m_b_zoom_reg_L1<-matrix(data = NA,nrow = nrow(mat1),ncol = ncol(mat1))
  m_b_zoom_reg_cor<-matrix(data = NA,nrow = nrow(mat1),ncol = ncol(mat1))
 
  
  dim(mat2)
  r1=3 
  # c<-1.5
  r2=7
  m=proc.time()
  
  
  for(i in (w1+1):(nrow(mat1)-w1))
  {
    for(j in (w1+1):(nrow(mat1)-w1))
    {
          # i<-40
          # j<-40
      N1<-c()
      for(k1 in (i-r1):(i+r1))
      {
        for(l1 in (j-r1):(j+r1))
        {
          if((k1-i)^2 + (l1-j)^2 <= (r1)^2)
          {
            N1<-c(N1,mat1[k1,l1])
          }
        }
      }
      msd<-c()
      mad<-c()
      r<-c()
      xz<-c()
      yz<-c()
      cov1<-0
      # cv<-c()
      sd1<-0
      # s1<-c()
      sd2<-0
      # s2<-c()
      
      for(m1 in (i-r2):(i+r2))
      {
        for(n1 in (j-r2):(j+r2))
        {
          if((m1-i)^2+(n1-j)^2 <= (r2)^2)
          {
            # m1<- (i-r2+5)
            # n1<- (j-r2+5)
            N2<-c()
            N<-0
            for(s1 in (m1-r1):(m1+r1))                              
            {
              for(t1 in (n1-r1):(n1+r1))
              {
                if((s1-m1)^2+(t1-n1)^2 <= (r1)^2)
                {
                  N2<-c(N2,mat2[s1,t1])
                  N<-N+1
                }
              }
            }
            cov1<-cov(N1,N2)
            # cv<-c(cv,cov1)
            
            sd1<-sqrt(var(N1))
            # s1<-c(s1,sd1)
            
            sd2<-sqrt(var(N2))
            # s2<-c(s2,sd2)
            r<-c(r,(cov1/(sd1*sd2)))
            msd<-c(msd,sum((N1-N2)^2)/N)
            mad<-c(mad,sum(abs(N1-N2))/N)
            xz<-c(xz,round(s1))
            yz<-c(yz,round(t1))
            
          }
        }
      }
      #r[is.infinite(r)]<-0
      index_L2<-which.min(msd)
      index_L1<-which.min(mad)
      index_cor<-which.max(r)
      
      m_b_zoom_reg_L2[round(xz[index_L2])-r1,round(yz[index_L2])-r1]<- mat1[i,j]
      m_b_zoom_reg_L1[round(xz[index_L1])-r1,round(yz[index_L1])-r1]<- mat1[i,j]
      m_b_zoom_reg_cor[round(xz[index_cor])-r1,round(yz[index_cor])-r1]<-mat1[i,j]
      
    }
  }
  
  n=proc.time()
  n-m
 typeof(setdiff(1:5,2:4))
  
  windows(10,10)
  par(mfrow=c(2,3))
  image(rot90(mat1[(w1+1):(nrow(m_b_zoom_reg_L1)-w1),(w1+1):(nrow(m_b_zoom_reg_L1)-w1)],3),col=gray(seq(0,1,length=256)))
  image(rot90(mat2[(w1+1):(nrow(m_b_zoom_reg_L1)-w1),(w1+1):(nrow(m_b_zoom_reg_L1)-w1)],3),col=gray(seq(0,1,length=256)))
  image(rot90(m_b_zoom_reg_L2[(w1+1):(nrow(m_b_zoom_reg_L1)-w1),(w1+1):(nrow(m_b_zoom_reg_L1)-w1)],3),col = grey(seq(0,1,length=256)))
  image(rot90(m_b_zoom_reg_L1[(w1+1):(nrow(m_b_zoom_reg_L1)-w1),(w1+1):(nrow(m_b_zoom_reg_L1)-w1)],3),col = grey(seq(0,1,length=256)))
  image(rot90(m_b_zoom_reg_cor[(w1+1):(nrow(m_b_zoom_reg_L2)-w1),(w1+1):(nrow(m_b_zoom_reg_L2)-w1)],3),col = grey(seq(0,1,length=256)))
  
  
  m_b_zoom_reg_L1[40,40]
  
  h<-3
  for(i in (w1+r1):(nrow(mat1)-w1-r1+1))
  {
    for(j in (w1+r1):(ncol(mat1)-w1-r1+1))
    {
      if(is.na(m_b_zoom_reg_L1[i,j]))
      {
        s<-c()
        for( k in (i-h):(i+h))
        {
          for(l in (j-h):(j+h))
          {
            s<- c(s,m_b_zoom_reg_L1[k,l]) 
          }
        }
        m_b_zoom_reg_L1[i,j]<-mean(s,na.rm = TRUE)
      }
      if(is.na(m_b_zoom_reg_L2[i,j]))
      {
        s1<-c()
        for( k in (i-h):(i+h))
        {
          for(l in (j-h):(j+h))
          {
            s1<- c(s1,m_b_zoom_reg_L2[k,l]) 
          }
        }
        m_b_zoom_reg_L2[i,j]<-mean(s1,na.rm = TRUE)
      }
      if(is.na(m_b_zoom_reg_cor[i,j]))
      {
        s2<-c()
        for( k in (i-h):(i+h))
        {
          for(l in (j-h):(j+h))
          {
            s2<- c(s2,m_b_zoom_reg_cor[k,l]) 
          }
        }
        m_b_zoom_reg_cor[i,j]<-mean(s2,na.rm = TRUE)
      }
    }
  }
  
  
  ##Finding Anomaly##
  
  anomaly_L1<-matrix(data=NA,nrow = nrow(zoom_pad),ncol = ncol(zoom_pad))
  anomaly_L2<-matrix(data=NA,nrow = nrow(zoom_pad),ncol = ncol(zoom_pad))
  anomaly_cor<-matrix(data=NA,nrow = nrow(zoom_pad),ncol = ncol(zoom_pad))
  
  for(i in (w1+r1):(nrow(zoom_pad)-w1-r1+1))
  {
    for(j in (w1+1+r1):(ncol(zoom_pad)-w1-r1+1))
    {
      if((abs(m_b_pad[i,j]-m_b_zoom_reg_L1[i,j])< (max(max(m_b_pad))-min(min(m_b_pad)))*0.15))
        anomaly_L1[i,j]<-0
      
      else
        anomaly_L1[i,j]<-1
      
      if((abs(m_b_pad[i,j]-m_b_zoom_reg_L2[i,j])< (max(max(m_b_pad))-min(min(m_b_pad)))*0.15))
        anomaly_L2[i,j]<-0
      
      else
        anomaly_L2[i,j]<-1
      
      if((abs(m_b_pad[i,j]-m_b_zoom_reg_cor[i,j])< (max(max(m_b_pad))-min(min(m_b_pad)))*0.15))
        anomaly_cor[i,j]<-0
      
      else
        anomaly_cor[i,j]<-1
      
      
      # if((zoom_pad[i,j]>0.5 && m_b_zoom_reg_L1[i,j]<0.5)|(zoom_pad[i,j]>0.5 && m_b_zoom_reg_L2[i,j]<0.5))
      #    #|(zoom_pad[i,j]>0.5 && m_b_zoom_reg_cor[i,j]<0.5)#)
      # {
      #   anomaly_L1[i,j]<-1
      #   anomaly_L2[i,j]<-1
      #   anomaly_cor[i,j]<-1
      # }
      # else if((zoom_pad[i,j]<0.5 && m_b_zoom_reg_L1[i,j]>0.5)|(zoom_pad[i,j]<0.5 && m_b_zoom_reg_L2[i,j]>0.5))
      #         #|(zoom_pad[i,j]<0.5 && m_b_zoom_reg_cor[i,j]>0.5))
      # {
      #   anomaly_L1[i,j]<-1
      #   anomaly_L2[i,j]<-1
      #   anomaly_cor[i,j]<-1
      # }
      # else
      # {
      #   anomaly_L1[i,j]<-0
      #   anomaly_L2[i,j]<-0
      #   anomaly_cor[i,j]<-0
      # } 
    }
  }
  
  windows(10,10)
  par(mfrow=c(1,3))
  image(anomaly_L1[(w1+1):(nrow(m_b_zoom_reg_L2)-w1),(w1+1):(nrow(m_b_zoom_reg_L2)-w1)],col = grey(seq(0,1,length=256)))
  image(anomaly_L2[(w1+1):(nrow(m_b_zoom_reg_L2)-w1),(w1+1):(nrow(m_b_zoom_reg_L2)-w1)],col = grey(seq(0,1,length=256)))
  image(anomaly_cor[(w1+1):(nrow(m_b_zoom_reg_L2)-w1),(w1+1):(nrow(m_b_zoom_reg_L2)-w1)],col = grey(seq(0,1,length=256)))
  
  anomaly[110:128,110:128]
  anomaly[110,111]
  
  anomaly_L1[(w1+1):(nrow(m_b_zoom_reg_L2)-w1),(w1+1):(nrow(m_b_zoom_reg_L2)-w1)]
  
  round(2.3)
  
  ##MSD between zoomed image and registered zoomed image##
  s1<-0
  s2<-0
  s3<-0
  n<-0
  for(i in (w1+r1):(nrow(mat2)-w1-r1+1))
  {
    for(j in (w1+r1):(nrow(mat2)-w1-r1+1))
    {
      
      s1<-s1+(mat2[i,j]-m_b_zoom_reg_L1[i,j])^2
      s2<-s2+(mat2[i,j]-m_b_zoom_reg_L2[i,j])^2
      s3<-s3+(mat2[i,j]-m_b_zoom_reg_cor[i,j])^2
      n<-n+1
    }
  }
  
  msd_L1<-sqrt(s1/n)
  msd_L2<-sqrt(s2/n)
  msd_cor<-sqrt(s3/n)
  
  c(msd_L1,msd_L2,msd_cor)
  
install.packages("magick")
library(magick)
test<- image_read("C:\\Users\\cssc\\OneDrive\\Documents\\cameraman_test_image.bmp")
test1<-image_read("C:\\Users\\cssc\\OneDrive\\Documents\\test_image_lena.png")
test

test[,,1]


install.packages("imager")
library(imager)

plot(test1)

test1[1:112,1:112]

remove(list = ls())

library(jpeg)

##Work on image 'lena'
test_lena<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\test_image_lena.jpg")

test_lena<-test_lena + matrix(rnorm(nrow(test_lena)*ncol(test_lena),0,0.01),nrow = nrow(test_lena),ncol = ncol(test_lena))


test_lena_low<-test_lena[32:307,32:345]

image(rot90(test_lena_low,k=3),col = grey(seq(0,1,length=256)))

test_lena_high<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\test_lena_higher_resolution.jpg")

test_lena_high<-test_lena_high + matrix(rnorm(nrow(test_lena)*ncol(test_lena),0,0.01),nrow = nrow(test_lena),ncol = ncol(test_lena))

library(matlab)
test_lena_pad<-padarray(test_lena,c(15,15),"symmetric","both")
test_lena_high_pad<-padarray(test_lena_high,c(15,15),"symmetric","both")

windows(10,10)
par(mfrow=c(1,2))
image(rot90(test_lena,k=3),col = grey(seq(0,1,length=256)))
image(rot90(test_lena_high,k=3),col = grey(seq(0,1,length=256)))



if(exists("rasterImage")){
  plot(1:2, type='n')
  rasterImage(test_lena,-1,-1,1,1)
}

str(test_lena_high)


##work on 'peppers' 256x256
test_pepper<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\test_images\\peppers_gray_jpg_256.jpg")
dim(test_pepper)
test_pepper_low<-test_pepper[64:191,64:191]

image(rot90(test_pepper_low,k=3),col = grey(seq(0,1,length=256)))

test_pepper_high<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\test_pepper_high.jpg")

windows(10,10)
par(mfrow=c(1,2))

image(rot90(test_pepper,k=3),col = grey(seq(0,1,length=256)))

image(rot90(test_pepper_high,k=3),col = grey(seq(0,1,length=256)))



#Work on image 'cameraman'

test_cam<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\cameraman_test_image.jpg")
test_cam<-test_cam + matrix(rnorm(nrow(test_cam)*ncol(test_cam),0,0.01),nrow = nrow(test_cam),ncol = ncol(test_cam))

dim(test_cam[1:145,56:200])

image(rot90(test_cam[1:145,56:200],k=3),col = grey(seq(0,1,length=256)))

test_cam_low<-test_cam[1:145,56:200]

test_cam_high<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\test_cam_high.jpg")
test_cam_high<-test_cam_high + matrix(rnorm(nrow(test_cam_high)*ncol(test_cam_high),0,0.01),nrow = nrow(test_cam_high),ncol = ncol(test_cam_high))


windows(10,10)
par(mfrow=c(1,2))
image(rot90(test_cam,k=3),col = grey(seq(0,1,length=256)))
image(rot90(test_cam_high,k=3),col = grey(seq(0,1,length=256)))



##Work on 'lena' original 512x512

lena_test<- readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\lena_original.jpg")
lena_test<- lena_test+ matrix(rnorm(nrow(lena_test)*ncol(lena_test),0,0.01),nrow = nrow(lena_test),ncol = ncol(lena_test))
dim(lena_test_low)
lena_test_low<-lena_test[64:448,64:448]
dim(lena_test_high)
image(rot90(lena_test_low,k=3),col=grey(seq(0,1,length=256)))
lena_test_high<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\lena_ori_high.jpg")
lena_test_high<-lena_test_high+ matrix(rnorm(nrow(lena_test)*ncol(lena_test),0,0.01),nrow = nrow(lena_test),ncol = ncol(lena_test))

par(mfrow=c(1,2))
image(rot90(lena_test,k=3),col=grey(seq(0,1,length=256)))
image(rot90(lena_test_high,k=3),col=grey(seq(0,1,length=256)))

#work on tile1
tile1<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\tile_1.jpg")
tile1_low<-tile1[17:112,17:112]
windows(10,10)
image(rot90(tile1_low,k=3),col=grey(seq(0,1,length=256)))
tile1_high<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\gray_mosaic_high.jpg")
dim(tile1)


#Work on tile2
tile2<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\tile_2.jpg")
tile2_low<-tile2[17:112,17:112]
windows(10,10)
image(rot90(tile2_high,k=3),col=grey(seq(0,1,length=256)))
tile2_high<-readJPEG("C:\\Users\\cssc\\OneDrive\\Documents\\tile2_high.jpg")
dim(tile2)

#Work on tile3
tile3<-readJPEG("tile4.jpg")
tile3_zm<-tile3[8:57,8:57]

jpeg(filename = "tile4_z.jpg",height=64,width=64)
par(bty="n",xaxt="n",yaxt="n",mar=c(0,0,0,0))
image(rot90(tile3_zm,3),col = grey(seq(0,1,length=256)))
dev.off()

windows(10,10)
image(rot90(tile3,k=3),col=grey(seq(0,1,length=256)))
tile3_z<-readJPEG("tile4_z.jpg")
dim(tile3_z)



## Registration algorithm for zoomed in real image##
w1<-20

img<- padarray(girl_zm,c(w1,w1),"symmetric","both")  #Zoomed image
img_zoom<-padarray(girl,c(w1,w1),"symmetric","both")     #Reference image
  
set.seed(8)
m=proc.time()
stop<-FALSE

dim(pepper_zm)

for(r2 in seq(2,5,by=1))
{
  
  for(r1 in c(r2-1,r2))
  {
    msd_L1<-NULL
    msd_L2<-NULL
    msd_cor<-NULL
    avg_msd_L1<-NULL
    avg_msd_L2<-NULL
    avg_msd_cor<-NULL
    var_msd_L1<-NULL
    var_msd_L2<-NULL
    var_msd_cor<-NULL
  for(x in 1:10)
  {
r1<-2
r2<-6

img<-img + matrix(rnorm(nrow(img)*ncol(img),0,0.01),nrow = nrow(img),ncol = ncol(img))
img_zoom<-img_zoom + matrix(rnorm(nrow(img)*ncol(img),0,0.01),nrow = nrow(img),ncol = ncol(img))

img_L2<-matrix(rep(NA,nrow(img)*ncol(img)),nrow = nrow(img),ncol = ncol(img))
img_L1<-matrix(rep(NA,nrow(img)*ncol(img)),nrow = nrow(img),ncol = ncol(img))
img_cor<-matrix(rep(NA,nrow(img)*ncol(img)),nrow = nrow(img),ncol = ncol(img))



#dim(test_lena_high_reg_cor)

for(i in (w1+1):(nrow(img)-w1))
{
  for(j in (w1+1):(ncol(img)-w1))
  {
    #i<-309
    #j=371
    
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
          #r<-c(r,(cov1/(sd1*sd2)))
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
    #r[is.infinite(r)]<-0
    index_L2<-which.min(msd)
    index_L1<-which.min(mad)
    index_cor<-which.max(r)
    img_L2[xz[index_L2]-r1,yz[index_L2]-r1]<- img[i,j]
    img_L1[xz[index_L1]-r1,yz[index_L1]-r1]<- img[i,j]
    img_cor[xz[index_cor]-r1,yz[index_cor]-r1]<- img[i,j]
  }
  if(stop){break}
}

# windows(10,10)
# par(mfrow=c(2,3))
# image(rot90(img,k=3),col=gray(seq(0,1,length=256)))
# image(rot90(img_zoom,k=3),col=gray(seq(0,1,length=256)))
# image(rot90(img_zoom_reg_L2,k=3),col = grey(seq(0,1,length=256)))
# image(rot90(img_zoom_reg_L1,k=3),col = grey(seq(0,1,length=256)))
# image(rot90(img_zoom_reg_cor,k=3),col = grey(seq(0,1,length=256)))

windows(10,10)
image(rot90(img_L1,k=3),col = grey(seq(0,1,length=256)))

# windows(10,10)
# par(mfrow=c(2,3))
# image(rot90(img[15:(nrow(img)-15),15:(ncol(img)-15)],k=3),col=gray(seq(0,1,length=256)))
# image(rot90(img_zoom[15:(nrow(img)-15),15:(ncol(img)-15)],k=3),col=gray(seq(0,1,length=256)))
# image(rot90(img_zoom_reg_L2[15:(nrow(img)-15),15:(ncol(img)-15)],k=3),col = grey(seq(0,1,length=256)))
# image(rot90(img_zoom_reg_L1[15:(nrow(img)-15),15:(ncol(img)-15)],k=3),col = grey(seq(0,1,length=256)))
# image(rot90(img_zoom_reg_cor[15:(nrow(img)-15),15:(ncol(img)-15)],k=3),col = grey(seq(0,1,length=256)))
# 
# dim(test_lena)


# max(max(test_lena))-min(min(test_lena))

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
n<-0
for(i in (w1+r1):(nrow(img)-w1-r1+1))
{
  for(j in (w1+r1):(nrow(img)-w1-r1+1))
  {
    if(is.na(img_L1[i,j])==FALSE)
    t1<-t1+(img_zoom[i,j]-img_L1[i,j])^2
    
    if(is.na(img_L2[i,j])==FALSE)
    t2<-t2+(img_zoom[i,j]-img_L2[i,j])^2
    
    if(is.na(img_cor[i,j])==FALSE)
    t3<-t3+(img_zoom[i,j]-img_cor[i,j])^2
    
    n<-n+1
  }
}


msd_L1<-c(msd_L1,sqrt(t1/n))
msd_L2<-c(msd_L2,sqrt(t2/n))
msd_cor<-c(msd_cor,sqrt(t3/n))


}
    assign(paste("avg_msd_L1",r1,r2,sep = "_"),mean(msd_L1))
    assign(paste("avg_msd_L2",r1,r2,sep = "_"),mean(msd_L2))
    assign(paste("avg_msd_cor",r1,r2,sep = "_"),mean(msd_cor))
    assign(paste("sd_msd_L1",r1,r2,sep = "_"),sqrt(var(msd_L1)))
    assign(paste("sd_msd_L2",r1,r2,sep = "_"),sqrt(var(msd_L2)))
    assign(paste("sd_msd_cor",r1,r2,sep = "_"),sqrt(var(msd_cor)))
    # avg_msd_L1<-c(avg_msd_L1,mean(msd_L1))
    # avg_msd_L2<-c(avg_msd_L2,mean(msd_L2))
    # avg_msd_cor<-c(avg_msd_cor,mean(msd_cor))
    # var_msd_L1<-c(var_msd_L1,var(msd_L1))
    # var_msd_L2<-c(var_msd_L2,var(msd_L2))
    # var_msd_cor<-c(var_msd_cor,var(msd_cor))
  }

}

n=proc.time()

m-n


for(i in 1:2)
{
  assign(paste("pg",i,sep = ""),i)
}

avg_msd_L1<-NULL
avg_msd_L2<-NULL
avg_msd_cor<-NULL
sd_msd_L1<-NULL
sd_msd_L2<-NULL
sd_msd_cor<-NULL

c<-1
while(c<=2)
{
  for(i in seq(1,3,by=0.2))
  {
    
    avg_msd_L1<-c(avg_msd_L1,print(get(paste("avg_msd_L1",i,i*c,sep = "_"))))
    avg_msd_L2<-c(avg_msd_L2,print(get(paste("avg_msd_L2",i,i*c,sep = "_"))))
    avg_msd_cor<-c(avg_msd_cor,print(get(paste("avg_msd_cor",i,i*c,sep = "_"))))
    sd_msd_L1<-c(sd_msd_L1,print(get(paste("sd_msd_L1",i,i*c,sep = "_"))))
    sd_msd_L2<-c(sd_msd_L2,print(get(paste("sd_msd_L2",i,i*c,sep = "_"))))
    sd_msd_cor<-c(sd_msd_cor,print(get(paste("sd_msd_cor",i,i*c,sep = "_"))))
  }
  c<-c+0.5
}

avg_msd_cor_1.2_1.2

matplot(cbind(avg_msd_L1,avg_msd_L2,avg_msd_cor),type = "l",lty = 1:3)

matplot(cbind(sd_msd_L1,sd_msd_L2,sd_msd_cor),type = "l",lty=1:3)


avg_sd<-data.frame(avg_msd_L1,avg_msd_L2,avg_msd_cor,sd_msd_L1,sd_msd_L2,sd_msd_cor)
View(avg_sd)

write.csv(avg_sd,file = "C:\\Users\\cssc\\OneDrive\\Documents\\avg_sd.csv")



img<- tile3
img_zoom<-tile3_high

img<-img + matrix(rnorm(nrow(img)*ncol(img),0,0.01),nrow = nrow(img),ncol = ncol(img))
img_zoom<-img_zoom + matrix(rnorm(nrow(img)*ncol(img),0,0.01),nrow = nrow(img),ncol = ncol(img))

img_zoom_reg_L2<-matrix(rep(0,nrow(img)*ncol(img)),nrow = nrow(img),ncol = ncol(img))
img_zoom_reg_L1<-matrix(rep(0,nrow(img)*ncol(img)),nrow = nrow(img),ncol = ncol(img))
img_zoom_reg_cor<-matrix(rep(0,nrow(img)*ncol(img)),nrow = nrow(img),ncol = ncol(img))
stop<-FALSE
r1<-2.8
c<-2
r2=c*r1

#dim(test_lena_high_reg_cor)


for(i in (r1+r2+1):(nrow(img)-r1-r2))
{
  for(j in (r1+r2+1):(ncol(img)-r1-r2))
  {
    #i<-309
    #j=371
    
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
    index_L2<-which.min(msd)
    index_L1<-which.min(mad)
    index_cor<-which.max(r)
    img_zoom_reg_L2[i,j]<- img_zoom[xz[index_L2],yz[index_L2]]
    img_zoom_reg_L1[i,j]<- img_zoom[xz[index_L1],yz[index_L1]]
    img_zoom_reg_cor[i,j]<- img_zoom[xz[index_cor],yz[index_cor]]
  }
  if(stop){break}
}
dim(img)


L1_img<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)
L1_img[(r1+1):((nrow(L1_img)-r1)),(r1+1):((ncol(L1_img)-r1))]<-img_L1[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

L2_img<--matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)
L2_img[(r1+1):((nrow(L2_img)-r1)),(r1+1):((ncol(L2_img)-r1))]<-img_L2[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

cor_img<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)
cor_img[(r1+1):((nrow(cor_img)-r1)),(r1+1):((ncol(cor_img)-r1))]<-img_cor[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]


windows(10,10)
jpeg(filename = "girl_L1(2,6).jpg",height=128,width=128)
par(bty="n",xaxt="n",yaxt="n",mai=c(0,0,0,0))
# image(rot90(img,k=3),col = grey(seq(0,1,length=256)))
#image(rot90(img_zoom,k=3),col = grey(seq(0,1,length=256)))
image(rot90(L1_img,k=3),col = grey(seq(0,1,length=256)))
box(which="figure")
# image(rot90(img_zoom_reg_L1,k=3),col = grey(seq(0,1,length=256)))
# image(rot90(img_zoom_reg_cor,k=3),col = grey(seq(0,1,length=256)))
dev.off()
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


L1_anmy<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)
L1_anmy[(r1+1):((nrow(L1_anmy)-r1)),(r1+1):((ncol(L1_anmy)-r1))]<-anomaly_L1[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

L2_anmy<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)
L2_anmy[(r1+1):((nrow(L2_anmy)-r1)),(r1+1):((ncol(L2_anmy)-r1))]<-anomaly_L2[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

cor_anmy<-matrix(NA,nrow = nrow(img)-2*w1,ncol = ncol(img)-2*w1)
cor_anmy[(r1+1):((nrow(cor_anmy)-r1)),(r1+1):((ncol(cor_anmy)-r1))]<-anomaly_cor[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]

windows(10,10)
jpeg(filename = "anomaly_pepper_cor(4,6).jpg",height=256,width=256)
par(bty="n",xaxt="n",yaxt="n",mai=c(0,0,0,0))
# par(mfrow=c(1,3))
# image(rot90(anomaly_L1,k=3),col = grey(seq(0,1,length=256)))
image(rot90(cor_anmy,k=3),col = grey(seq(0,1,length=256)))
# image(rot90(anomaly_cor,k=3),col = grey(seq(0,1,length=256)))
dev.off()
0.6*qnorm(0.9)


length(c(img_cor[(w1+r1+1):((nrow(img)-w1-r1)),(w1+r1+1):((ncol(img)-w1-r1))]))

t=0
for(i in (w1+r1+1):((nrow(img)-w1-r1)))
{
  for(j in (w1+r1+1):((ncol(img)-w1-r1)))
  {
    if(is.na(img_cor[i,j])==FALSE)
      t=t+1
  }
}


