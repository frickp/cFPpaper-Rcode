basefile	<-	'~/Dropbox/Shared Vito/cFP paper/Figures/'
#Code to automatically add PC compatibility
basefile2	<-	gsub('\\\\','/',path.expand(basefile))
xx			<-	try(setwd(basefile2))
if(class(xx)=="try-error"){basefile2 <- paste('C:/Users/Peter',substring(basefile,2),sep="")}
read.dir 	<- paste(basefile2,"Data used for figures (pulled by R)",sep="")
write.dir	<- paste(basefile2,"Figure parts",sep="")
setwd(read.dir)

##################################################################
# Pull all data
##################################################################

library('fracprolif')
source('cFP norm (72h).r')

# subset data to use only CHX500 and control
# exclude any data with a non-zero value for FSK, TRM, or X17A

myIDlist	<-	c(	"B04_J_E0_CHX500_TRM0_X17A0_Plate1",
					"B03_A_E0_CHX500_TRM0_X17A0_Plate1",
					"B05_H_E0_CHX500_TRM0_X17A0_Plate1"
				)

cfp.ctrl	<-	cfp[apply((apply(cfp[,c('Erl','CHX','FSK','TRM','X17A')], c(1,2), FUN=function(x){x==0})),1,FUN=all),]
ctrl.rate 	<-	coef(lm(nl2 ~ Time.day, data=cfp.ctrl))['Time.day']/24	# doubling rate in hours

cfp			<-	subset(cfp, id %in% myIDlist)
cfp1		<-	subset(cfp, id==myIDlist[1])


# split_tot =
#function (t, x0, q0, d1, d2, q1, q2, ad1, aq1, ad2, aq2, br = 20) 
#{
#    t = time
#    x0 = starting number in prolif compartment
#    q0 = starting number in quiesc compartment
#    d1 = initial div rate
#    d2 = second div rate
#    q1 = initial qui rate
#    q2 = second qui rate
#    ad1 = initial death (apoptosis) rate from div compartment
#    aq1 = initial death (apoptosis) rate from qui compartment
#    ad2 = second death (apoptosis) rate from div compartment
#    aq2 = second death (apoptosis) rate from qui compartment
#    br = time to switch between rates (break)
#
#    usually keep ad1==aq1 and ad2==aq2
#    
#    qbr <- prolif_q(br, x0, q0, d1, q1, ad1, aq1)
#    dbr <- prolif_d(br, x0, q0, d1, q1, ad1, aq1)
#    before <- prolif_tot(t, x0, q0, d1, q1, ad1, aq1)
#    after <- prolif_tot(t - br, dbr, qbr, d2, q2, ad2, aq2)
#    answer <- before
#    answer[t > br] <- after[t > br]
#    answer
#}

x	<-	cfp1$Time.day

m1.fit		<-	nls(	
	algorithm = "port",
	cfp1$l2 ~ split_tot(	x, 
								cfp1$l2[1],   # x0=fixed
								0,            # no cells in qui compartment at start
								ctrl.rate,    # fixed
								div.rate2,    # floating
								0,            # fixed
								qr2, # floating
								0,            # ad1=fixed
								0,            # aq1=fixed
								death.rate,   # floating
								death.rate,   # floating
								br
				) ,          # floating
		start=list(					div.rate2=ctrl.rate/8, 
									qr2=0.01,
									death.rate=0.0001,
									br=2
		),
		lower=list(					div.rate2=0.0001, 
									qr2=0.0001,
									death.rate=0.00001,
									br=0
		),
		upper=list(					div.rate2=ctrl.rate, 
									qr2=0.04,
									death.rate=0.00015,
									br=4
		),
		trace=TRUE,
		data=cfp1
)



tot		<-	split_tot(	t=x,
						x0=cfp1$l2[1], 
						q0=0, 
						d1=ctrl.rate,
						d2=coef(m1.fit)['div.rate2'], 
						q1=0,
						q2=coef(m1.fit)['qr2'], 
						ad1=0, 
						aq1=0, 
						ad2=coef(m1.fit)['death.rate'], 
						aq2=coef(m1.fit)['death.rate'], 
						br=coef(m1.fit)['br'])


plot(l2 ~ Time.day, data=cfp1)
lines(x,tot)
						
						
				

