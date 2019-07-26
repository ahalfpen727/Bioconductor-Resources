#code is based on this animation https://gist.github.com/jebyrnes/b34930da0052a86f5ffe254ce9900357
#and this visualisation on this data http://www.realclimate.org/images/nsidc4.txt

#1. load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(animation)

#2. get data in right format
col.names <- c("Year", "Month", "Day", "Extent", "Missing", "Source_Data")

#3. Get Data: I downloaded the data from these ftps
#ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/NH_seaice_extent_final_v2.csv
#ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/NH_seaice_extent_nrt_v2.csv

nisdc1 = read.csv("NH_seaice_extent_final_v2.csv",skip=2,col.names=col.names)
nisdc2 = read.csv("NH_seaice_extent_nrt_v2.csv",skip=2,col.names=col.names)

nisdc <- data.frame(Year=c(nisdc1$Year,nisdc2$Year),
                    Month=c(nisdc1$Month,nisdc2$Month),
                    Day=c(nisdc1$Day,nisdc2$Day),
                    Extent=c(nisdc1$Extent,nisdc2$Extent),
                    Missing=c(nisdc1$Missing,nisdc2$Missing))



#get a day from each month. Early in dataset they only record every second day so this is a kludge
stripped <- filter(nisdc, Day == 1)
stripped4 <- filter(nisdc, Year==1988 & Month==1 & Day ==13) #month 1 of 1989 is wonky
stripped2 <- filter(nisdc, Day == 2) #get the months without a day 1
stripped3 <- filter(stripped2, Year<1988) #before 1989 the only every odd day issue happens
stripped5 <- rbind(stripped, stripped3)
stripped5 <- rbind(stripped5, stripped4)

mo <- months(seq(as.Date("1910/1/1"), as.Date("1911/1/1"), "months"))
mo <- gsub("(^...).*", "\\1", mo)


#4. make animation
saveGIF({

  for(i in 1979:2016){
 print(ggplot(stripped5 %>% filter(Year <= i),
           aes(x=Month, y=Extent, color=Year, group=Year)) +
        geom_line()+
        scale_color_gradient(low="blue", high="black", limits=c(1978, 2016), guide="none")+
        annotate(x=2, y=2.5, geom="text", label=i, size = 6) +
        annotate(x=1, y=17, geom="text", label = "10^6*km^2",parse = TRUE, size = 4) +
                ggtitle("Arctic Sea Ice Size 1978-2016")+
        scale_x_continuous(labels=mo, breaks=1:13) +
        scale_y_continuous(breaks=0:16) +
         theme_classic()
         +
        annotate(x=10, y=2, geom="text", label = "Chart by @iamreddave | Data source: nsidc.org", size = 3)

         )}
}, interval=0.3,ani.width = 600, ani.height = 600)
