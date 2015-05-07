
/* -*- Mode: C; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
/*
 *  Unit SGP4SDP4
 *           Author:  Dr TS Kelso 
 * Original Version:  1991 Oct 30
 * Current Revision:  1992 Sep 03
 *          Version:  1.50 
 *        Copyright:  1991-1992, All Rights Reserved 
 *
 *   Ported to C by:  Neoklis Kyriazis  April 10  2001
 *   Reentrancy mods by Alexandru Csete OZ9AEC
 *
 * 	 Ported to JavaScript by: Igor Kim 31.03.2015
*/

/*
*	 Модель SGP4/SDP4 Автор Dr TS Kelso
*	 Реализация математической модели SGP4 на языке JavaScript
*	 Портировано с языка C (библиотека программы GPredict): Игорь Ким 31.03.2015
*/
SGDP4 = function(orbital_elements)
{
	//Анализ орбитальных элементов из TLE
	var epy = Number(orbital_elements["epoch_year"])
	//Пересчет года в полный формат
	if(epy<57)
	{
		this.epoch_year = epy+2000
	}
	else
	{
		this.epoch_year = epy+1900
	}

	this.bstar_mantissa = Number((orbital_elements["bstar_mantissa"]*1e-5).toFixed(5))
	this.bstar_exponent = Number("1e" + orbital_elements["bstar_exponent"])
	this.bstar  		= this.bstar_mantissa*this.bstar_exponent

	this.epoch 			= orbital_elements["epoch"]
	this.xincl 			= orbital_elements["inclination"]*torad
	this.xnodeo 		= orbital_elements["right_ascension"]*torad
	this.eo 			= Number((orbital_elements["eccentricity"]*1e-7).toFixed(7))
	this.omegao 	 	= orbital_elements["argument_of_perigee"]*torad
	this.xmo 			= orbital_elements["mean_anomaly"]*torad
	this.xno 			= orbital_elements["mean_motion"]*twopi/1440.0

	this.a1 			= Math.pow(xke/this.xno,tothrd)
	this.cosio 			= Math.cos(this.xincl)
	this.theta2 		= this.cosio*this.cosio
	this.x3thm1 		= 3*this.theta2-1.0
	this.eosq 			= this.eo*this.eo
	this.betao2 		= 1-this.eosq
	this.betao 			= Math.sqrt(this.betao2)
	this.del1 			= 1.5*ck2*this.x3thm1/(this.a1*this.a1*this.betao*this.betao2)
	this.ao 			= this.a1*(1-this.del1*((1.0/3.0)+this.del1*(1.0+(134.0/81.0)*this.del1)))
	this.delo 			= 1.5*ck2*this.x3thm1/(this.ao*this.ao*this.betao*this.betao2)
	this.xnodp 			= this.xno/(1.0+this.delo) //original_mean_motion
	this.aodp 			= this.ao/(1.0-this.delo) //semi_major_axis

	this.orbital_period	= twopi / this.xnodp

	//deep space check
	this.DEEP_SPACE_FLAG = (this.orbital_period >= 225.0)

	// For perigee less than 220 kilometers, the "simple" flag is set
	// and the equations are truncated to linear variation in sqrt a
	// and quadratic variation in mean anomaly.  Also, the c3 term,
	// the delta omega term, and the delta m term are dropped.
	this.SIMPLE_SGP_FLAG = ((this.aodp*(1.0-this.eo)/ae) < (220.0/xkmper+ae))&&(!this.DEEP_SPACE_FLAG)

	// For perigee below 156 km, the values of s and qoms2t are altered.
	var s4 		= s
	var qoms24 	= qoms2t
	this.perige = (this.aodp*(1.0-this.eo)-ae)*xkmper
	this.apoge 	= (this.aodp*(1.0+this.eo)-ae)*xkmper
	if (this.perige < 156.0) {
		s4 		= (this.perige <= 98.0)?20:this.perige-78.0
		qoms24	= Math.pow(((120.0-s4)*ae/xkmper),4)

		s4 		= s4/xkmper+ae
	}

	if (this.DEEP_SPACE_FLAG) {
		this.sing 		= Math.sin(this.omegao)
		this.cosg 		= Math.cos(this.omegao)
	}
	var pinvsq 		= 1.0/(this.aodp*this.aodp*this.betao2*this.betao2)
	var tsi 		= 1.0/(this.aodp-s4)
	this.eta 		= this.aodp*this.eo*tsi
	var etasq 		= this.eta*this.eta
	var eeta 		= this.eo*this.eta
	var psisq 		= Math.abs(1.0-etasq)
	var coef 		= qoms24*Math.pow(tsi,4)
	var coef1 		= coef/Math.pow(psisq,3.5)
	var c2 			= coef1*this.xnodp*(this.aodp*(1.0+1.5*etasq+eeta*(4.0+etasq))+
						0.75*ck2*tsi/psisq*this.x3thm1*(8.0+3.*etasq*(8.0+etasq)))
	this.c1 		= this.bstar*c2
	this.sinio 		= Math.sin(this.xincl)
	var a3ovk2 		= -xj3/ck2*Math.pow(ae,3)
	this.x1mth2 	= 1.0-this.theta2
	this.c4 		= 2.0*this.xnodp*coef1*this.aodp*this.betao2*(this.eta*(2.0+0.5*etasq)+
						this.eo*(0.5+2.0*etasq)-2.0*ck2*tsi/(this.aodp*psisq)*
						(-3.0*this.x3thm1*(1.0-2.0*eeta+etasq*(1.5-0.5*eeta))+
							0.75*this.x1mth2*(2.0*etasq-eeta*(1.0+etasq))*Math.cos((2.0*this.omegao)))
						)
	this.c5 		= 2.0*coef1*this.aodp*this.betao2*(1.0+2.75*(etasq+eeta)+eeta*etasq)
	this.theta4 	= this.theta2*this.theta2
	var temp1 		= 3.0*ck2*pinvsq*this.xnodp
	var temp2 		= temp1*ck2*pinvsq
	var temp3 		= 1.25*ck4*pinvsq*pinvsq*this.xnodp
	this.xmdot 		= this.xnodp+0.5*temp1*this.betao*this.x3thm1+0.0625*temp2*this.betao*(13.0-78.0*this.theta2+137.0*this.theta4)
	var x1m5th 		= 1.0-5.0*this.theta2
	this.omgdot 	= -0.5*temp1*x1m5th+0.0625*temp2*(7.0-114.0*this.theta2+395.0*this.theta4)+
						temp3*(3.0-36.0*this.theta2+49.0*this.theta4)
	var xhdot1 		= -temp1*this.cosio
	this.xnodot 	= xhdot1+(0.5*temp2*(4.0-19.0*this.theta2)+2.0*temp3*(3.0-7.0*this.theta2))*this.cosio
	this.xnodcf 	= 3.5*this.betao2*xhdot1*this.c1
	this.t2cof 		= 1.5*this.c1
	this.xlcof 		= 0.125*a3ovk2*this.sinio*(3.0+5.0*this.cosio)/(1.0+this.cosio)
	this.aycof		= 0.25*a3ovk2*this.sinio
	this.x7thm1 	= 7.0*this.theta2-1.0

	if (this.SIMPLE_SGP_FLAG != 1) {
		this.sinmo 	= Math.sin(this.xmo)
		this.xmcof 	= -tothrd*coef*this.bstar*ae/eeta
		this.delmo 	= Math.pow((1.0+this.eta*Math.cos(this.xmo)),3)

		var c3 		= coef*tsi*a3ovk2*this.xnodp*ae*this.sinio/this.eo //??
		this.omgcof = this.bstar*c3*Math.cos(this.omegao)
		var c1sq 	= this.c1*this.c1
		this.d2 	= 4.*this.aodp*tsi*c1sq
		var temp 	= this.d2*tsi*this.c1/3.0
		this.d3 	= (17.0*this.aodp+s4)*temp
		this.d4 	= 0.5*temp*this.aodp*tsi*(221.0*this.aodp+31.*s4)*this.c1
		this.t3cof 	= this.d2+2.0*c1sq
		this.t4cof 	= 0.25*(3.0*this.d3+this.c1*(12.0*this.d2+10.0*c1sq))
		this.t5cof 	= 0.2*(3.0*this.d4+12.0*this.c1*this.d3+6.0*this.d2*this.d2+15.*c1sq*(2.0*this.d2+c1sq))
	}
	if (this.DEEP_SPACE_FLAG) {
		this.deep('dpinit')
	}
}

//Вычислительная функция
SGDP4.prototype.calc = function(ts, tsince) {
	var tsince = (typeof tsince !== 'undefined')?tsince:elapsedTime1(this.epoch_year,this.epoch,ts)
	//UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
	var xmdf 	= this.xmo+this.xmdot*tsince
	this.omgadf = this.omegao+this.omgdot*tsince
	var xnoddf 	= this.xnodeo+this.xnodot*tsince
	var tsq 	= tsince*tsince
	this.xnode 	= xnoddf+this.xnodcf*tsq
	var tempa 	= 1.0-this.c1*tsince
	var tempe 	= this.bstar*this.c4*tsince
	var templ	= this.t2cof*tsq

	if (this.DEEP_SPACE_FLAG) {
		this.xn 	= this.xnodp
		// Update for deep-space secular effects
		this.xll 	= xmdf
		this.t 		= tsince
		this.deep('dpsec')
		xmdf 		= this.xll
		a 			= Math.pow(xke / this.xn, tothrd) * tempa * tempa
		this.em 	= this.em - tempe
		xmam 		= xmdf + this.xnodp * templ

		// Update for deep-space periodic effects
		this.xll 	= xmam
		this.deep('dpper')
		xmam 		= this.xll
		xl 			= xmam + this.omgadf +this.xnode
	} else {
		var xmp 	= xmdf		
		if (this.SIMPLE_SGP_FLAG != 1) {
			var delomg 	= this.omgcof*tsince
			var delm 	= this.xmcof*(Math.pow((1.0+this.eta*Math.cos(xmdf)),3)-this.delmo)
			var temp 	= delomg+delm
			var xmp 	= xmdf+temp
			this.omgadf = omgadf-temp
			var tcube 	= tsq*tsince
			var tfour 	= tsince*tcube
			var tempa 	= tempa-this.d2*tsq-this.d3*tcube-this.d4*tfour
			var tempe 	= tempe+this.bstar*this.c5*(Math.sin(xmp)-this.sinmo)
			var templ 	= templ+this.t3cof*tcube+tfour*(this.t4cof+tsince*this.t5cof)
		}
		var a 		= this.aodp*tempa*tempa
		this.em 	= this.eo-tempe
		var xl 		= xmp+this.omgadf+this.xnode+this.xnodp*templ		
	}
	var beta 	= Math.sqrt(1.0-this.em*this.em)
	this.xn 	= xke/Math.pow(a,1.5)

	// long period periodics
	var axn 	= this.em*Math.cos(this.omgadf)
	var temp 	= 1.0/(a*beta*beta)
	var xll 	= temp*this.xlcof*axn
	var aynl 	= temp*this.aycof
	var xlt 	= xl+xll
	var ayn 	= this.em*Math.sin(this.omgadf)+aynl
	//console.log('xlt: '+xl+' '+xll+' = '+xlt)
	// solve keplers equation
	var capu 	= fmod2p(xlt-this.xnode)
	//console.log('capu: '+xlt+' '+this.xnode+' = '+capu)
	var temp2 	= capu
	
	i = 0
	do {
		sinepw = Math.sin(temp2)
		cosepw = Math.cos(temp2)
		temp3 = axn * sinepw
		temp4 = ayn * cosepw
		temp5 = axn * cosepw
		temp6 = ayn * sinepw
		epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2

		//console.log(i+': '+temp2+' '+sinepw+', '+cosepw)
		if (Math.abs(epw - temp2) <= e6a) {
			break
		}
		temp2 = epw
	} while (i++ < 10)

	/*
	for (i=1; i<=10; i++)
	{
		var sinepw 	= Math.sin(temp2)
		var cosepw 	= Math.cos(temp2)
		var temp3 	= axn*sinepw
		var temp4 	= ayn*cosepw
		var temp5 	= axn*cosepw
		var temp6 	= ayn*sinepw
		var epw 	= (capu-temp4+temp3-temp2)/(1.0-temp5-temp6)+temp2
		if (Math.abs(epw-temp2) <= e6a)
		{
			break
		}
		temp2 		= epw
	}
	*/

	//short period preliminary quantities
	var ecose 		= temp5+temp6
	var esine 		= temp3-temp4
	var elsq 		= axn*axn+ayn*ayn
	var temp 		= 1.0-elsq
	var pl 			= a*temp
	var r 			= a*(1.0-ecose)
	var temp1 		= 1.0/r
	var rdot 		= xke*Math.sqrt(a)*esine*temp1
	var rfdot 		= xke*Math.sqrt(pl)*temp1
	var temp2 		= a*temp1
	var betal 		= Math.sqrt(temp)
	var temp3 		= 1.0/(1.0+betal)
	var cosu 		= temp2*(cosepw-axn+ayn*esine*temp3)
	var sinu 		= temp2*(sinepw-ayn-axn*esine*temp3)

	//console.log(temp2+' '+cosepw+' '+' '+axn+' '+ayn+' '+esine+' '+temp3)
	//console.log(cosu+' '+sinu)
	/*TODO:
	 Check this equations 	
	 In GPredict: u = AcTan (sinu, cosu);
	*/
	/*var u 			= Math.atan2(sinu,cosu)
	if (u<0)
	{
		u+= 2* Math.PI
	}*/
	/****/
	u = Math.atan2(sinu, cosu)


	var sin2u 		= 2.0*sinu*cosu
	var cos2u 		= 2.0*cosu*cosu-1.0
	var temp 		= 1.0/pl
	var temp1 		= ck2*temp
	var temp2 		= temp1*temp

	// update for short periodics
	var rk 			= r*(1.0-1.5*temp2*betal*this.x3thm1)+0.5*temp1*this.x1mth2*cos2u
	var uk 			= u-0.25*temp2*this.x7thm1*sin2u
	var xnodek 		= this.xnode+1.5*temp2*this.cosio*sin2u
	var xinck 		= this.xinc+1.5*temp2*this.cosio*this.sinio*cos2u
	var rdotk 		= rdot-this.xn*temp1*this.x1mth2*sin2u
	var rfdotk 		= rfdot+this.xn*temp1*(this.x1mth2*cos2u+1.5*this.x3thm1)

	//console.log("rk = "+' '+r+' '+temp2+' '+betal+' '+temp1+' '+cos2u+'  = '+rk)
	//console.log('xinck = '+xinck)
	//orientation vectors
	var sinuk 		= Math.sin(uk)
	var cosuk 		= Math.cos(uk)
	var sinik 		= Math.sin(xinck)
	var cosik 		= Math.cos(xinck)
	var sinnok 		= Math.sin(xnodek)
	var cosnok 		= Math.cos(xnodek)
	var xmx 		= -sinnok*cosik
	var xmy 		= cosnok*cosik
	var ux 			= xmx*sinuk+cosnok*cosuk
	var uy 			= xmy*sinuk+sinnok*cosuk
	var uz 			= sinik*sinuk
	var vx 			= xmx*cosuk-cosnok*sinuk
	var vy 			= xmy*cosuk-sinnok*sinuk
	var vz 			= sinik*cosuk

	//console.log("ux = "+xmx+'*'+sinuk+'+'+cosnok+'*'+cosuk+'='+ux)
	//console.log("uy = "+xmy+'*'+sinuk+'+'+sinnok+'*'+cosuk+'='+uy)
	//console.log("uz = "+sinik+'*'+sinuk+'='+uz)

	// Position and velocity
	this.x 			= rk*ux
	this.y 			= rk*uy
	this.z 			= rk*uz
	this.xdot 		= rdotk*ux+rfdotk*vx
	this.ydot 		= rdotk*uy+rfdotk*vy
	this.zdot 		= rdotk*uz+rfdotk*vz
	// *******END OF SGP4

	/*
	Gpredict also returned 
	sat->phase = xlt - xnode - omgadf + twopi;
	if (sat->phase < 0)
		sat->phase += twopi;
	sat->phase = FMod2p (sat->phase);

	sat->tle.omegao1 = omega;
	sat->tle.xincl1  = xinck;
	sat->tle.xnodeo1 = xnodek;
	*/
	var result = {
		x	: this.x,
		y 	: this.y,
		z 	: this.z,
		xdot: this.xdot,
		ydot: this.ydot,
		zdot: this.zdot
	}
	return result
}

// DEEP
// This function is used by SDP4 to add lunar and solar
// perturbation effects to deep-space orbit objects.
SGDP4.prototype.deep = function(mode)
{
	// Entrance for deep space initialization
	if (mode === 'dpinit') {
		thetag_temp = thetag1(this.epoch,this.epoch_year)
		this.thgr 	= thetag_temp.THETAG
		this.ds50 	= thetag_temp.DS50

		eq 			= this.eo
		this.xnq 	= this.xnodp
		aqnv 		= 1.0/this.aodp
		this.xqncl 	= this.xincl
		xmao 		= this.xmo
		xpidot 		= this.omgdot + this.xnodot
		sinq 		= Math.sin(this.xnodeo)
		cosq 		= Math.cos(this.xnodeo)
		this.omegaq = this.omegao
		this.preep 	= 0

		// Initialize lunar solar terms
		day 		= this.ds50 + 18261.5 //Days since 1900 Jan 0.5
		if (day != this.preep) {
			preep 		= day
			xnodce 		= 4.5236020 - 9.2422029e-4 * day
			stem 		= Math.sin(xnodce)
			ctem 		= Math.cos(xnodce)
			this.zcosil = 0.91375164 - 0.03568096 * ctem
			this.zsinil = Math.sqrt(1.0 - this.zcosil * this.zcosil)
			this.zsinhl = 0.089683511 * stem / this.zsinil
			this.zcoshl = Math.sqrt(1.0 - this.zsinhl * this.zsinhl)
			c 			= 4.7199672 + 0.22997150 * day
			gam 		= 5.8351514 + 0.0019443680 * day
			this.zmol 	= fmod2p(c - gam)
			zx 			= 0.39785416 * stem / this.zsinil
			zy 			= this.zcoshl * ctem + 0.91744867 * this.zsinhl * stem
			zx 			= Math.atan2(zx,zy)
			zx 			= gam + zx - xnodce
			this.zcosgl = Math.cos(zx)
			this.zsingl = Math.sin(zx)
			this.zmos 	= 6.2565837 + 0.017201977 * day
			this.zmos 	= fmod2p(this.zmos)
		}// End if(day != preep)

		// Do solar terms
		this.savtsn = 1e20
		zcosg 		= zcosgs
		zsing 		= zsings
		zcosi 		= zcosis
		zsini 		= zsinis
		zcosh 		= cosq
		zsinh 		= sinq
		cc 			= c1ss
		zn 			= zns
		ze 			= zes
		zmo 		= this.zmos
		xnoi 		= 1.0 / this.xnq

		this.LUNAR_TERMS_DONE_FLAG = false

		// Loop breaks when Solar terms are done a second
		// time, after Lunar terms are initialized
		for(;;)  {
			// Solar terms done again after Lunar terms are done
			a1 	= zcosg * zcosh + zsing * zcosi * zsinh
			a3 	= -zsing * zcosh + zcosg * zcosi * zsinh
			a7 	= -zcosg * zsinh + zsing * zcosi * zcosh
			a8 	= zsing * zsini
			a9 	= zsing * zsinh + zcosg * zcosi * zcosh
			a10 = zcosg * zsini
			a2 	= this.cosio * a7 + this.sinio * a8
			a4 	= this.cosio * a9 + this.sinio * a10
			a5 	= -this.sinio * a7 + this.cosio * a8
			a6 	= -this.sinio*a9+ this.cosio*a10
			x1 	= a1*this.cosg+a2*this.sing
			x2 	= a3*this.cosg+a4*this.sing
			x3 	= -a1*this.sing+a2*this.cosg
			x4 	= -a3*this.sing+a4*this.cosg
			x5 	= a5*this.sing
			x6 	= a6*this.sing
			x7 	= a5*this.cosg
			x8 	= a6*this.cosg
			z31 = 12*x1*x1-3*x3*x3
			z32 = 24*x1*x2-6*x3*x4
			z33 = 12*x2*x2-3*x4*x4
			z1 	= 3*(a1*a1+a2*a2)+z31*this.eosq
			z2 	= 6*(a1*a3+a2*a4)+z32*this.eosq
			z3 	= 3*(a3*a3+a4*a4)+z33*this.eosq
			z11 = -6*a1*a5+this.eosq*(-24*x1*x7-6*x3*x5)
			z12 = -6*(a1*a6+a3*a5)+ this.eosq*(-24*(x2*x7+x1*x8)-6*(x3*x6+x4*x5))
			z13 = -6*a3*a6+this.eosq*(-24*x2*x8-6*x4*x6)
			z21 = 6*a2*a5+this.eosq*(24*x1*x5-6*x3*x7)
			z22 = 6*(a4*a5+a2*a6)+ this.eosq*(24*(x2*x5+x1*x6)-6*(x4*x7+x3*x8))
			z23 = 6*a4*a6+this.eosq*(24*x2*x6-6*x4*x8)
			z1 	= z1+z1+this.betao2*z31
			z2 	= z2+z2+this.betao2*z32
			z3 	= z3+z3+this.betao2*z33
			s3 	= cc*xnoi
			s2 	= -0.5*s3/this.betao
			s4 	= s3*this.betao
			s1 	= -15*eq*s4
			s5 	= x1*x3+x2*x4
			s6 	= x2*x3+x1*x4
			s7 	= x2*x4-x1*x3
			se 	= s1*zn*s5
			si 	= s2*zn*(z11+z13)
			sl 	= -zn*s3*(z1+z3-14-6*this.eosq)
			sgh = s4*zn*(z31+z33-6)
			sh 	= -zn*s2*(z21+z23)
			if (this.xqncl < 5.2359877e-2) {
				sh = 0;
			}

			this.ee2 	= 2*s1*s6
			this.e3 	= 2*s1*s7
			this.xi2 	= 2*s2*z12
			this.xi3 	= 2*s2*(z13-z11)
			this.xl2 	= -2*s3*z2
			this.xl3 	= -2*s3*(z3-z1)
			this.xl4 	= -2*s3*(-21-9*this.eosq)*ze
			this.xgh2 	= 2*s4*z32
			this.xgh3 	= 2*s4*(z33-z31)
			this.xgh4 	= -18*s4*ze
			this.xh2 	= -2*s2*z22
			this.xh3 	= -2*s2*(z23-z21)

			if (this.LUNAR_TERMS_DONE_FLAG) {
				break
			}

			// Do lunar terms
			this.sse 	= se
			this.ssi 	= si
			this.ssl 	= sl
			this.ssh 	= sh/this.sinio
			this.ssg 	= sgh-this.cosio*this.ssh
			this.se2 	= this.ee2
			this.si2 	= this.xi2
			this.sl2 	= this.xl2
			this.sgh2 	= this.xgh2
			this.sh2 	= this.xh2
			this.se3 	= this.e3
			this.si3 	= this.xi3
			this.sl3 	= this.xl3
			this.sgh3 	= this.xgh3
			this.sh3 	= this.xh3
			this.sl4 	= this.xl4
			this.sgh4 	= this.xgh4
			zcosg 		= this.zcosgl
			zsing 		= this.zsingl
			zcosi 		= this.zcosil
			zsini 		= this.zsinil
			zcosh 		= this.zcoshl*cosq+this.zsinhl*sinq
			zsinh 		= sinq*this.zcoshl-cosq*this.zsinhl
			zn 			= znl
			cc 			= c1l
			ze 			= zel
			zmo 		= this.zmol
			this.LUNAR_TERMS_DONE_FLAG = true
		} // End of for(;;)

		this.sse 	= this.sse+se
		this.ssi 	= this.ssi+si
		this.ssl 	= this.ssl+sl
		this.ssg 	= this.ssg+sgh-this.cosio/this.sinio*sh
		this.ssh 	= this.ssh+sh/this.sinio

		// Geopotential resonance initialization for 12 hour orbits */
		//TODO: true or false??
		//sat->flags &= ~RESONANCE_FLAG;
		//sat->flags &= ~SYNCHRONOUS_FLAG;
		this.RESONANCE_FLAG 	= false
		this.SYNCHRONOUS_FLAG 	= false

		if( !((this.xnq < 0.0052359877) && (this.xnq > 0.0034906585)) )  {
			if( (this.xnq < 0.00826) || (this.xnq > 0.00924) ) {
				return
			}

			if (eq < 0.5) {
				return
			}
			this.RESONANCE_FLAG = true
			eoc 				= eq*this.eosq
			g201 				= -0.306-(eq-0.64)*0.440
			if (eq <= 0.65)  {
				g211 = 3.616-13.247*eq+16.290*this.eosq
				g310 = -19.302+117.390*eq-228.419* this.eosq+156.591*eoc
				g322 = -18.9068+109.7927*eq-214.6334* this.eosq+146.5816*eoc
				g410 = -41.122+242.694*eq-471.094* this.eosq+313.953*eoc
				g422 = -146.407+841.880*eq-1629.014* this.eosq+1083.435*eoc
				g520 = -532.114+3017.977*eq-5740* this.eosq+3708.276*eoc
			} else  {
				g211 = -72.099+331.819*eq-508.738*this.eosq+266.724*eoc
				g310 = -346.844+1582.851*eq-2415.925*this.eosq+1246.113*eoc
				g322 = -342.585+1554.908*eq-2366.899*this.eosq+1215.972*eoc
				g410 = -1052.797+4758.686*eq-7193.992*this.eosq+3651.957*eoc
				g422 = -3581.69+16178.11*eq-24462.77*this.eosq+ 12422.52*eoc
				if (eq <= 0.715) {
					g520 = 1464.74-4664.75*eq+3763.64*this.eosq
				} else {
					g520 = -5149.66+29936.92*eq-54087.36*this.eosq+31324.56*eoc
				}
			} /* End if (eq <= 0.65) */

			if (eq < 0.7) {
				g533 = -919.2277+4988.61*eq-9064.77*this.eosq+5542.21*eoc
				g521 = -822.71072+4568.6173*eq-8491.4146*this.eosq+5337.524*eoc
				g532 = -853.666+4690.25*eq-8624.77*this.eosq+ 5341.4*eoc
			} else {
				g533 = -37995.78+161616.52*eq-229838.2*this.eosq+109377.94*eoc
				g521 = -51752.104+218913.95*eq-309468.16*this.eosq+146349.42*eoc
				g532 = -40023.88+170470.89*eq-242699.48*this.eosq+115605.82*eoc
			} /* End if (eq <= 0.7) */

			sini2 		= this.sinio*this.sinio
			f220 		= 0.75*(1+2*this.cosio+this.theta2)
			f221 		= 1.5*sini2
			f321 		= 1.875*this.sinio*(1-2*this.cosio-3*this.theta2)
			f322 		= -1.875*this.sinio*(1+2*this.cosio-3*this.theta2)
			f441 		= 35*sini2*f220
			f442 		= 39.3750*sini2*sini2
			f522 		= 9.84375*this.sinio*(sini2*(1-2*this.cosio-5*this.theta2)+0.33333333*(-2+4*this.cosio+6*this.theta2))
			f523 		= this.sinio*(4.92187512*sini2*(-2-4*this.cosio+10*this.theta2)+6.56250012*(1+2*this.cosio-3*this.theta2))
			f542 		= 29.53125*this.sinio*(2-8*this.cosio+this.theta2*(-12+8*this.cosio+10*this.theta2))
			f543 		= 29.53125*this.sinio*(-2-8*this.cosio+this.theta2*(12+8*this.cosio-10* this.theta2))
			xno2 		= this.xnq*this.xnq
			ainv2 		= aqnv*aqnv
			temp1 		= 3*xno2*ainv2
			temp 		= temp1*root22
			this.d2201 	= temp*f220*g201
			this.d2211 	= temp*f221*g211
			temp1 		= temp1*aqnv
			temp 		= temp1*root32
			this.d3210 	= temp*f321*g310
			this.d3222 	= temp*f322*g322
			temp1 		= temp1*aqnv
			temp 		= 2*temp1*root44
			this.d4410 	= temp*f441*g410
			this.d4422 	= temp*f442*g422
			temp1 		= temp1*aqnv
			temp 		= temp1*root52
			this.d5220 	= temp*f522*g520
			this.d5232 	= temp*f523*g532
			temp 		= 2*temp1*root54
			this.d5421 	= temp*f542*g521
			this.d5433 	= temp*f543*g533
			this.xlamo 	= xmao+this.xnodeo+this.xnodeo-this.thgr-this.thgr
			bfact 		= this.xmdot+this.xnodot+this.xnodot-thdt-thdt
			bfact 		= bfact+this.ssl+this.ssh+this.ssh
		} /* if( !(this.xnq < 0.0052359877) && (this.xnq > 0.0034906585) ) */
		else {
			this.RESONANCE_FLAG 	= true
			this.SYNCHRONOUS_FLAG 	= true
			/* Synchronous resonance terms initialization */
			g200 		= 1+this.eosq*(-2.5+0.8125*this.eosq)
			g310 		= 1+2*this.eosq
			g300 		= 1+this.eosq*(-6+6.60937*this.eosq)
			f220 		= 0.75*(1+this.cosio)*(1+this.cosio)
			f311 		= 0.9375*this.sinio*this.sinio*(1+3*this.cosio)-0.75*(1+this.cosio)
			f330 		= 1+this.cosio
			f330 		= 1.875*f330*f330*f330
			this.del1 	= 3*this.xnq*this.xnq*aqnv*aqnv
			this.del2 	= 2*this.del1*f220*g200*q22
			this.del3 	= 3*this.del1*f330*g300*q33*aqnv
			this.del1 	= this.del1*f311*g310*q31*aqnv
			this.fasx2 	= 0.13130908
			this.fasx4 	= 2.8843198
			this.fasx6 	= 0.37448087
			this.xlamo 	= xmao+this.xnodeo+this.omegao-this.thgr
			bfact 		= this.xmdot+xpidot-thdt
			bfact 		= bfact+this.ssl+this.ssg+this.ssh
		} /* End if( !(xnq < 0.0052359877) && (xnq > 0.0034906585) ) */

		this.xfact 	= bfact-this.xnq

		/* Initialize integrator */
		this.xli 	= this.xlamo
		this.xni 	= this.xnq
		this.atime 	= 0
		this.stepp 	= 720
		this.stepn 	= -720
		this.step2 	= 259200
		/* End case dpinit: */
		return;
	}
	/* Entrance for deep space secular effects */
	else if (mode === 'dpsec')
	{
		this.xll 		= this.xll+this.ssl*this.t
		this.omgadf 	= this.omgadf+this.ssg*this.t
		this.xnode 		= this.xnode+this.ssh*this.t
		this.em 		= this.eo+this.sse*this.t
		this.xinc 		= this.xincl+this.ssi*this.t
		if (this.xinc < 0) {
			this.xinc = -this.xinc
			this.xnode = this.xnode + pi
			this.omgadf = this.omgadf-pi
		}
		//if( ~sat->flags & RESONANCE_FLAG ) return
		if (!this.RESONANCE_FLAG) {
			return false
		}

		do {
			if ( (this.atime == 0) || ((this.t >= 0) && 
				(this.atime < 0)) || ((this.t < 0) && (this.atime >= 0)) ) {
				/* Epoch restart */
				delt = ( this.t >= 0 )?this.stepp:this.stepn

				this.atime 	= 0
				this.xni 	= this.xnq
				this.xli 	= this.xlamo
			} else {	  
				if( Math.abs(this.t) >= Math.abs(this.atime) ) {
					delt = ( this.t > 0 )?this.stepp:this.stepn
				}
			}

			do {
				if ( Math.abs(this.t-this.atime) >= this.stepp ) {
					this.DO_LOOP_FLAG 		= true
					this.EPOCH_RESTART_FLAG	= false
					/*sat->flags |= DO_LOOP_FLAG;
					sat->flags &= ~EPOCH_RESTART_FLAG;*/
				}
				else {
					ft 						= this.t-this.atime
					//sat->flags &= ~DO_LOOP_FLAG;
					this.DO_LOOP_FLAG 		= false
				}

				if( Math.abs(this.t) < Math.abs(this.atime) ) {
					if (this.t >= 0) {
						delt = this.stepn
					} else {
						delt = this.stepp
					}
					//sat->flags |= (DO_LOOP_FLAG | EPOCH_RESTART_FLAG);
					this.DO_LOOP_FLAG 		= true
					this.EPOCH_RESTART_FLAG	= true					
				}

				/* Dot terms calculated */
				//if (sat->flags & SYNCHRONOUS_FLAG) 
				if (this.SYNCHRONOUS_FLAG = true) {
					xndot = this.del1*Math.sin(this.xli-this.fasx2)+this.del2*Math.sin(2*(this.xli-this.fasx4))
								+this.del3*Math.sin(3*(this.xli-this.fasx6))
					xnddt = this.del1*Math.cos(this.xli-this.fasx2)+2*this.del2*Math.cos(2*(this.xli-this.fasx4))
								+3*this.del3*Math.cos(3*(this.xli-this.fasx6))
				} else {
					xomi 	= this.omegaq+this.omgdot*this.atime
					x2omi 	= xomi+xomi
					x2li 	= this.xli+this.xli
					xndot 	= this.d2201*Math.sin(x2omi+this.xli-g22)+this.d2211*Math.sin(this.xli-g22)
								+this.d3210*Math.sin(xomi+this.xli-g32)+this.d3222*Math.sin(-xomi+this.xli-g32)
								+this.d4410*Math.sin(x2omi+x2li-g44)+this.d4422*Math.sin(x2li-g44)+this.d5220*Math.sin(xomi+this.xli-g52)
								+this.d5232*Math.sin(-xomi+this.xli-g52)+this.d5421*Math.sin(xomi+x2li-g54)
								+this.d5433*Math.sin(-xomi+x2li-g54)
					xnddt = this.d2201*Math.cos(x2omi+this.xli-g22)+this.d2211*Math.cos(this.xli-g22)
								+this.d3210*Math.cos(xomi+this.xli-g32)+this.d3222*Math.cos(-xomi+this.xli-g32)
								+this.d5220*Math.cos(xomi+this.xli-g52)+this.d5232*Math.cos(-xomi+this.xli-g52)
								+2*(this.d4410*Math.cos(x2omi+x2li-g44)+this.d4422*Math.cos(x2li-g44)
								+this.d5421*Math.cos(xomi+x2li-g54)+this.d5433*Math.cos(-xomi+x2li-g54))

				} /* End of if (isFlagSet(SYNCHRONOUS_FLAG)) */
				xldot = this.xni+this.xfact
				xnddt = xnddt*xldot

				//if(sat->flags & DO_LOOP_FLAG) 
				if (this.DO_LOOP_FLAG) {
					this.xli = this.xli+xldot*delt+xndot*this.step2
					this.xni = this.xni+xndot*delt+xnddt*this.step2
					this.atime = this.atime+delt
				}
			}
			while ((this.DO_LOOP_FLAG) &&(!this.EPOCH_RESTART_FLAG))
		}
		while ((this.DO_LOOP_FLAG) && (this.EPOCH_RESTART_FLAG))

		this.xn = this.xni+xndot*ft+xnddt*ft*ft*0.5
		xl 		= this.xli+xldot*ft+xndot*ft*ft*0.5
		temp 	= -this.xnode+this.thgr+this.t*thdt

		//if (~sat->flags & SYNCHRONOUS_FLAG)
		if (!this.SYNCHRONOUS_FLAG) {
			this.xll = xl+temp+temp
		} else {
			this.xll = xl-this.omgadf+temp
		}

		return
		/*End case dpsec: */
	}
	/* Entrance for lunar-solar periodics */
	else if (mode === 'dpper') {
		sinis = Math.sin(this.xinc)
		cosis = Math.cos(this.xinc)
		if (Math.abs(this.savtsn-this.t) >= 30) {
			this.savtsn 	= this.t
			zm 				= this.zmos+zns*this.t
			zf 				= zm+2*zes*Math.sin(zm)
			sinzf 			= Math.sin(zf)
			f2 				= 0.5*sinzf*sinzf-0.25
			f3 				= -0.5*sinzf*Math.cos(zf)
			ses 			= this.se2*f2+this.se3*f3
			sis 			= this.si2*f2+this.si3*f3
			sls 			= this.sl2*f2+this.sl3*f3+this.sl4*sinzf
			this.sghs 		= this.sgh2*f2+this.sgh3*f3+this.sgh4*sinzf
			this.shs 		= this.sh2*f2+this.sh3*f3
			zm 				= this.zmol+znl*this.t
			zf 				= zm+2*zel*Math.sin(zm)
			sinzf 			= Math.sin(zf)
			f2 				= 0.5*sinzf*sinzf-0.25
			f3 				= -0.5*sinzf*Math.cos(zf)
			sel 			= this.ee2*f2+this.e3*f3
			sil 			= this.xi2*f2+this.xi3*f3
			sll 			= this.xl2*f2+this.xl3*f3+this.xl4*sinzf
			this.sghl 		= this.xgh2*f2+this.xgh3*f3+this.xgh4*sinzf
			this.sh1 		= this.xh2*f2+this.xh3*f3
			this.pe 		= ses+sel
			this.pinc 		= sis+sil
			this.pl 		= sls+sll
		}
		pgh 		= this.sghs+this.sghl
		ph 			= this.shs+this.sh1
		this.xinc 	= this.xinc+this.pinc
		this.em 	= this.em+this.pe

		if (this.xqncl >= 0.2) {
			/* Apply periodics directly */
			ph 			= ph/this.sinio
			pgh 		= pgh-this.cosio*ph
			this.omgadf = this.omgadf+pgh
			this.xnode 	= this.xnode+ph
			this.xll 	= this.xll+this.pl
		} else {
			/* Apply periodics with Lyddane modification */
			sinok 		= Math.sin(this.xnode)
			cosok 		= Math.cos(this.xnode)
			alfdp 		= sinis*sinok
			betdp		= sinis*cosok
			dalf 		= ph*cosok+this.pinc*cosis*sinok
			dbet 		= -ph*sinok+this.pinc*cosis*cosok
			alfdp 		= alfdp+dalf
			betdp 		= betdp+dbet
			this.xnode 	= fmod2p(this.xnode)
			xls 		= this.xll+this.omgadf+cosis*this.xnode
			dls 		= this.pl+pgh-this.pinc*this.xnode*sinis
			xls 		= xls+dls
			xnoh 		= this.xnode
			this.xnode 	= Math.atan2(alfdp,betdp)

			/* This is a patch to Lyddane modification */
			/* suggested by Rob Matson. */
			if(Math.abs(xnoh-this.xnode) > pi) {
				if(this.xnode < xnoh) {
					this.xnode +=twopi
				} else {
					this.xnode -=twopi
				}
			}

			this.xll 		= this.xll+this.pl
			this.omgadf 	= xls-this.xll-Math.cos(this.xinc)*this.xnode
		} /* End case dpper: */
		return
	} else {
		return 0
	}
}