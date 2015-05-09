/*
	Дополнительные функции и константы
*/
//Константы WGS84
var WGS84_a 	= 6378137 //Semi-major axis of the Earth in meters
var WGS84_b 	= 6356752 //Semi-minor axis of the Earth in meters
var WGS84_1of 	= 298.26 //Flattening Factor of the Earth (1/f)
var WGS84_f 	= 1/WGS84_1of //Earth's flattening term in WGS-84 (= 1/298.26)
var WGS84_e 	= Math.sqrt(2*WGS84_f - WGS84_f*WGS84_f)

//Константы, необходимые для вычисления
var torad 	= Math.PI/180 //Перевод в радианы
var pi 		= 3.14159265 //π
var pio2 	= 1.57079633 //π/2
var x3pio2 	= 4.71238898 //3π/2
var twopi 	= 6.2831853 //2π
var e6a 	= 1.0e-6 //10^-6
var tothrd 	= 0.66666667 //2/3
var xj2 	= 1.082616e-3
var xj3 	= -2.53881e-6
var xj4 	= -1.65597e-6
var xke 	= 0.743669161e-1
var xkmper	= WGS84_a/1000 // Earth's equational radius in WGS-84 (km)
var minpday = 1440.0 // Minutes per day
var xmnpda 	= minpday
var ae 		= 1.0 // distance units/Earth radii
var ck2 	= 5.413080e-4
var ck4 	= 6.2098875e-7
var qo 		= ae +120.0/xkmper
var qoms2t 	= 1.88027916e-9
var s 		= 1.0+78.0/xkmper
var de2ra 	= 0.174532925e-1 //radians/degree .174532925E-1
var eflat 	= WGS84_f //Earth's flattening term in WGS-84 (= 1/298.26)

/*GPredict const*/
var secday 	=  8.6400e4        /* Seconds per day */
var omega_E =  1.0027379
var omega_ER= 6.3003879
var zcosis  = 9.1744867e-1
var zsinis  = 3.9785416e-1
var zsings  = -9.8088458e-1
var zcosgs  = 1.945905e-1
var zcoshs  = 1
var zsinhs  = 0
var zns     = 1.19459e-5
var c1ss    = 2.9864797e-6
var zes     = 1.675E-2
var znl     = 1.5835218e-4
var c1l     = 4.7968065e-7
var zel     = 5.490E-2
var q22     = 1.7891679e-6
var q31     = 2.1460748e-6
var q33     = 2.2123015e-7
var g22     = 5.7686396
var g32     = 9.5240898e-1
var g44     = 1.8014998
var g52     = 1.0508330
var g54     = 4.4108898
var root22  = 1.7891679e-6
var root32  = 3.7393792e-7
var root44  = 7.3636953e-9
var root52  = 1.1428639e-7
var root54  = 2.1765803e-9
var thdt    = 4.3752691e-3
var rho     = 1.5696615e-1
var mfactor = 7.292115e-5


function dateFromTimestamp(ts) {
	dt = new Date()
	dt.setTime(ts)

	return dt
}
function sortZ(a,b) {
	if (a.z < b.z)
		return -1
	if (a.z > b.z)
		return 1
	return 0
}
function xdiff(points) {
	for(var index in points) { 
		if (points.hasOwnProperty(index)) {
			if (typeof max === 'undefined') max = points[index].x
			if (typeof min === 'undefined') min = points[index].x
			
			if (points[index].x > max) max = points[index].x
			if (points[index].x < min) min = points[index].x
		}
	}
	return (typeof max !== 'undefined' && typeof min !== 'undefined')?Math.abs(max-min):0
}
function ydiff(points) {
	for(var index in points) { 
		if (points.hasOwnProperty(index)) {
			if (typeof max === 'undefined') max = points[index].y
			if (typeof min === 'undefined') min = points[index].y
			
			if (points[index].y > max) max = points[index].y
			if (points[index].y < min) min = points[index].y
		}
	}
	return (typeof max !== 'undefined' && typeof min !== 'undefined')?Math.abs(max-min):0	
}
function majorMinorAxis(x,y) {
	return {semiMajor:x, semiMinor:y}
}

function footprint_calc(satpos) {
	satpos.footprint = 2.0 * xkmper * Math.acos(xkmper/satpos.w)
	return satpos
}

function coverage_calc (satpos, mode) {
    warped 	= false;
    numrc 	= 1;
    /* Range circle calculations.
     * Borrowed from gsat 0.9.0 by Xavier Crehueras, EB3CZS
     * who borrowed from John Magliacane, KD2BD.
     * Optimized by Alexandru Csete and William J Beksi.
     */
	ssplat = satpos.lat * torad
	ssplon = satpos.lng * torad
	beta = (0.5 * satpos.footprint) / xkmper
	points 			= []
	mirrorPoints 	= []
	for (azi = 0; azi <= 180; azi++) {
		azimuth = azi*torad
		rangelat = Math.asin(Math.sin(ssplat) * Math.cos(beta) + Math.cos(azimuth) * Math.sin(beta) * Math.cos(ssplat))
		num = Math.cos(beta) - (Math.sin(ssplat) * Math.sin(rangelat))
		dem = Math.cos(ssplat) * Math.cos(rangelat)
		    
		if (azi == 0 && (beta > pio2 - ssplat)) {
			rangelon = ssplon + pi
		}
		    
		else if (azi == 180 && (beta > pio2 + ssplat)) {
			rangelon = ssplon + pi
		}
		        
		else if (Math.abs(num / dem) > 1.0) {
			rangelon = ssplon
		}
		        
		else {
			if ((180 - azi) >= 0) {
				rangelon = ssplon - acos2(num, dem)
			} else {
				rangelon = ssplon + acos2(num, dem)
			}
		}

		while (rangelon < -pi) {
			rangelon += twopi
		}

		while (rangelon > (pi)) {
			rangelon -= twopi
		}

		rangelat = rangelat / de2ra
		rangelon = rangelon / de2ra
		mlon = mirror_lon(satpos, rangelon)

		if (azi==0) {
			left1 	= Cesium.Cartesian3.fromDegrees(rangelon, rangelat, 0)
		} 
		if (azi==90) {
			left2 		= Cesium.Cartesian3.fromDegrees(rangelon, rangelat, 0)
			right2 		= Cesium.Cartesian3.fromDegrees(mlon, rangelat, 0)
			semiMinor 	= Cesium.Cartesian3.distance(left2,right2)
		}

		if (azi==180) {
			right1 		= 	Cesium.Cartesian3.fromDegrees(mlon, rangelat, 0)	
			semiMajor 	= Cesium.Cartesian3.distance(left1,right1)
		}

		if (typeof mode !== 'undefined') {
			points.push(rangelon, rangelat)
			mirrorPoints.push(mlon, rangelat)			
		} else {
			point = {lng: rangelon, lat: rangelat}
			point = Cesium.Cartesian3.fromDegrees(rangelon, rangelat, 0)	
			points.push(point)	

			point = {lng: mlon, lat: rangelat}
			point = Cesium.Cartesian3.fromDegrees(mlon, rangelat, 0)
			points.push(point)
		}
    }
    if (typeof mode !== 'undefined') {
    	if (semiMinor > semiMajor) {
			return {semiMajor: semiMinor, semiMinor: semiMajor}
    	} else {
    		return {semiMajor: semiMajor, semiMinor: semiMinor}
    	}
	} else {
		points.sort(sortZ)
	}
    return points
}

function getFlagForKeyCode(keyCode) {
	switch (keyCode) {
	case 'Q'.charCodeAt(0):
		return 'moveForward';
	case 'E'.charCodeAt(0):
		return 'moveBackward';
	case 'W'.charCodeAt(0):
		return 'moveUp';
	case 'S'.charCodeAt(0):
		return 'moveDown';
	case 'D'.charCodeAt(0):
		return 'moveRight';
	case 'A'.charCodeAt(0):
		return 'moveLeft';
	case 'T'.charCodeAt(0):
		return 'rotateLeft';
	case 'U'.charCodeAt(0):
		return 'rotateRight';
	case 'Y'.charCodeAt(0):
		return 'rotateLeft';
	case 'H'.charCodeAt(0):
		return 'rotateRight';
	default:
		return undefined;
	}
}

function TLEParse(line1, line2) {
	var orbital_elements= {
		//Элементы первой строки
		line_number_1 					: Number(line1.slice(0,0)),
		catalog_no_1 					: line1.slice(2,7),
		security_classification 		: Number(line1.slice(8,8)),
		international_identification 	: String(line1.slice(9,17)),
		epoch_year 						: Number(line1.slice(18,20)),
		epoch 							: Number(line1.substring(20,32)),
		first_derivative_mean_motion 	: Number(line1.substring(33,43)),
		second_derivative_mean_motion 	: Number(line1.substring(44,52)),
		bstar_mantissa					: Number(line1.substring(53,59)),
		bstar_exponent 					: Number(line1.substring(59,61)),
		ephemeris_type 					: Number(line1.substring(62,63)),
		element_number 					: Number(line1.substring(64,68)),
		check_sum_1 					: Number(line1.substring(69,69)),
		//Элементы второй строки
		line_number_2 					: Number(line2.slice(0,0)),
		catalog_no_2 					: Number(line2.slice(2,7)),
		inclination 					: Number(line2.substring(8,16)),
		right_ascension 				: Number(line2.substring(17,25)),
		eccentricity 					: Number(line2.substring(26,33)),
		argument_of_perigee 			: Number(line2.substring(34,42)),
		mean_anomaly 					: Number(line2.substring(43,51)),
		mean_motion 					: Number(line2.substring(52,63)),
		rev_number_at_epoch 			: Number(line2.substring(64,68)),
		check_sum_2 					: Number(line2.substring(68,69))
	}
	return orbital_elements
}


//Ephemeris JS Translation JuliadDay from ts
function JulianDay_calc(ts) {
	d = new Date();
	d.setTime(ts)

	var year 		= d.getUTCFullYear()
	var month 		= d.getUTCMonth()+1
	var day 		= d.getUTCDate()
	var hour 		= d.getUTCHours()
	var minute 		= d.getUTCMinutes()
	var second 		= d.getUTCSeconds()
	
	var calender 	= ""

	if(month <= 2)
	{
		var year 	= year - 1;
		var month 	= month + 12;
	} 

	var julian_day = Math.floor(365.25*(year+4716))+Math.floor(30.6001*(month+1))+day-1524.5
						+ hour/24 + minute/(24*60) + second/(24*60*60)

	if (calender == "julian")
	{ 
		var transition_offset = 0
	}
	else if (calender == "gregorian")
	{
		var tmp 				= Math.floor(year/100)
		var transition_offset 	= 2-tmp+Math.floor(tmp/4)
	}
	else if (julian_day<2299160.5)
	{
		var transition_offset = 0
	}
	else
	{
		var tmp 				= Math.floor(year/100)
		var transition_offset 	= 2-tmp+Math.floor(tmp/4)
	}

	var jd 	= julian_day+transition_offset

	return jd
}

//Ephemeris JS Translation GMST convertion from timestamp
function GMST_calc(ts){
	//load default values
	d = new Date()
	d.setTime(ts)
	var hours		= d.getUTCHours()
	var minutes		= d.getUTCMinutes()
	var seconds		= d.getUTCSeconds()
	var time_in_sec = hours*3600+minutes*60+seconds
	var jd 			= JulianDay_calc(ts)
	var rad 		= Math.PI/180

	//gmst at 0:00
	var t 				= (jd-2451545.0)/36525
	var gmst_at_zero 	= (24110.5484 + 8640184.812866*t+0.093104*t*t+0.0000062*t*t*t)/3600
	if(gmst_at_zero>24)
	{
		gmst_at_zero 	= gmst_at_zero%24
	}
	this.gmst_at_zero 	= gmst_at_zero

	//gmst at target time
	var gmst 			= gmst_at_zero+(time_in_sec * 1.00273790925)/3600

	//mean obliquity of the ecliptic
	e 					= 23+26.0/60 + 21.448/3600 -46.8150/3600*t -0.00059/3600*t*t +0.001813/3600*t*t*t
	//nutation in longitude
	omega 				= 125.04452-1934.136261*t+0.0020708*t*t+t*t*t/450000
	long1 				= 280.4665 + 36000.7698*t
	long2 				= 218.3165 + 481267.8813*t
	phai 				= -17.20*Math.sin(omega*rad)-(-1.32*Math.sin(2*long1*rad))-
								0.23*Math.sin(2*long2*rad) + 0.21*Math.sin(2*omega*rad)
	gmst 				= gmst + ((phai/15)*(Math.cos(e*rad)))/3600

	if(gmst<0) {
		gmst 	= gmst%24+24
	}

	if(gmst>24) {
		gmst 	= gmst%24
	}
	return gmst
}

//Сортировка элемента select
function sortSelect(selector) {
	var wrapper = document.getElementById(selector),
	nodes = wrapper.getElementsByTagName("option"),
	len = nodes.length,
	sorted = [];
	while (nodes[0]) 
	{
		sorted.push(new String(nodes[0].value));
		sorted[sorted.length-1].element = nodes[0];
		wrapper.removeChild(nodes[0]);
	}
	sorted = sorted.sort();
	for (var i = 0; i < len; i++) {
		wrapper.appendChild(sorted[i].element);
	}
}

//Функция возвращает текущее время в формате timestamp
function curTimeStampUTC() {
	date 				= new Date()
	return Date.UTC(date.getUTCFullYear(), date.getUTCMonth(), date.getUTCDate(), date.getUTCHours(), date.getUTCMinutes(), date.getUTCSeconds(), date.getUTCMilliseconds())	
}

//Асинхронная загрузка содержимого файла
function loadTextFileAjaxSync(filePath, mimeType, sendData) {
	var xmlhttp=new XMLHttpRequest()
	xmlhttp.open("POST",filePath,false)
	if (mimeType != null) {
		if (xmlhttp.overrideMimeType) 
		{
			xmlhttp.overrideMimeType(mimeType)
		}
	}
	xmlhttp.send(sendData)
	if (xmlhttp.status==200)
	{
		return xmlhttp.responseText
	}
	else {
		// TODO Throw exception
		return null
	}
}

//Создает HTMK-элемент element с текстом text
function createElement(element, text) {
	var newElement 	= document.createElement(element)
	newElement.text = text
	return newElement
}


//Отображение миллисекунд в формате 00:00.0
function formatSeconds(msec) {
	tsec 	= msec

	sec 	= Math.floor(tsec/1000)
	if (sec > 0) {
		msec = tsec%1000
	} else	{
		sec = 0
	}

	min 	= Math.floor(sec/60)
	if (min > 0) {
		sec = sec%60
	} else {
		min = 0
	}

	hour 	= Math.floor(min/60)
	if (hour > 0) {
		min = min%60
	} else {
		hour = 0
	}
	return hour+':'+('0'+min).slice(-2)+':'+('0'+sec).slice(-2)+'.'+msec
}

/* The function Julian_Date_of_Year calculates the Julian Date  */
/* of Day 0.0 of {year}. This function is used to calculate the */
/* Julian Date of any date by using Julian_Date_of_Year, DOY,   */
/* and Fraction_of_Day. */
function Julian_Date_of_Year(year) {	
	/* Astronomical Formulae for Calculators, Jean Meeus, */
	/* pages 23-25. Calculate Julian Date of 0.0 Jan year */
	year 	= year - 1
	i 		= parseInt(year / 100)
	A 		= i
	i 		= parseInt(A / 4)
	B 		= parseInt(2 - A + i)
	i 		= parseInt(365.25 * year)
	i 		+= parseInt(30.6001 * 14)
	jdoy 	= i + 1720994.5 + B
	return jdoy;
}

function getJulianDay_SatEpoch(year, dSatelliteEpoch) {
	year = year % 100
	year += (year < 57)?2000:1900

	dResult 	= Julian_Date_of_Year(year)
	dResult 	+= dSatelliteEpoch

	return dResult
}

function elapsedTime1(epoch_year,epoch,date) {
	satJD 		= getJulianDay_SatEpoch(epoch_year, epoch)
	curJD 		= JulianDay_calc(date)

	return (curJD-satJD)*minpday
}

//Время которое прошло с времени date
function elapsedTime(epoch_year,epoch,date) {
	var d = new Date()
	d.setTime(date)

	var year 			= d.getUTCFullYear()
	var month 			= d.getUTCMonth()+1
	var day				= d.getUTCDate()
	var hours 			= d.getUTCHours()
	var minutes 		= d.getUTCMinutes()
	var seconds 		= d.getUTCSeconds()
	var year2 			= epoch_year-1

	//Количество секунд до указанной даты в date
	var now_sec 		= Date.UTC(year, month-1, day, hours, minutes, seconds)
	//Количество секунд до дня запуска
	var epoch_sec 		= Date.UTC(year2, 11, 31, 0, 0, 0)+(epoch*24*60*60*1000)
	//Разница между запуском и указанной датой в переменной date
	var elapsed_time	= (now_sec-epoch_sec)/(60*1000)

	return elapsed_time.toFixed(2)
}

//Округление до определенного количества десятичных знаков
function round(arg, decimal_places) {
	dec 	= Math.pow(10, decimal_places)
	return Math.round(arg * dec) / dec
}

function thetag1(day, year) {
	/* Reference:  The 1992 Astronomical Almanac, page B6. */
	// double year,day,UT,jd,TU,GMST,_ThetaG;

	/* Modification to support Y2K */
	/* Valid 1957 through 2056     */
	//UT = fmod($day, $day);
	var jd = Julian_Date_of_Year(year) + day
	var TU = (jd - 2451545.0) / 36525
	var GMST = 24110.54841 + TU * (8640184.812866 + TU * (0.093104 - TU * 6.2e-6))
	GMST = GMST%secday
	var DS50 = jd - 2433281.5

	var result 	= {
		THETAG 	: fmod2p(6.3003880987 * DS50 + 1.72944494),
		DS50 	: DS50
	}
	return result
}

/* The function ThetaG calculates the Greenwich Mean Sidereal Time */
/* for an epoch specified in the format used in the NORAD two-line */
/* element sets. It has now been adapted for dates beyond the year */
/* 1999, as described above. The function ThetaG_JD provides the   */
/* same calculation except that it is based on an input in the     */
/* form of a Julian Date. */
function thetag(EP) {
	var YR = (EP + 2.0e-7) * 1.0e-3
	var JY = parseInt(YR)

	YR = JY

	var D = EP - YR * 1.0e3
	if(JY < 10) {
		JY += 80
	}

	var N = (JY - 69) / 4
	if(JY < 70) {
		var N = (JY - 72) / 4
	}

	var DS50 	= 7305.0 + 365.0 * (JY-70) + N + D
	var THETA 	= 1.72944494 + 6.3003880987 * DS50
	var TEMP 	= THETA / twopi
	var I 	= parseInt(TEMP)

	TEMP = I

	var THETAG = THETA - TEMP * twopi
	if(THETAG < 0.0) {
		THETAG += twopi
	}

	var result = {
		THETAG 	: THETAG,
		DS50 	: DS50
	}
	return result
}

/* Returns fractional part of double argument */
function Frac(arg) {
	return arg - Math.floor(arg);
}

/** \brief Arccosine implementation. 
 *
 * Returns a value between zero and two pi.
 * Borrowed from gsat 0.9 by Xavier Crehueras, EB3CZS.
 * Optimized by Alexandru Csete.
 */
function acos2(x, y) {
	if (x && y) {
		if (y > 0.0) {
			return Math.acos(x/y)
		}
		else if (y < 0.0) {
			return pi + Math.acos(x/y)
		}
	}
    return 0.0
}

/** \brief Mirror the footprint longitude. */
function mirror_lon(satpos, rangelon) {
	mlon = 0
	if (satpos.lng < 0.0) {
		/* western longitude */
		if (rangelon < 0.0) {
			/* rangelon has not been warped over */
			mlon = satpos.lng + Math.abs(rangelon - satpos.lng)
		}
		else {
			/* rangelon has been warped over */
			diff 	= 360.0 + satpos.lng - rangelon
			mlon 	= satpos.lng + diff
		}
	} else {
        /* eastern longitude */
		mlon = satpos.lng + Math.abs(rangelon - satpos.lng)
		if (mlon > 180.0) {
			mlon  -= 360
		}
    }
    return mlon
}

function thetag_JD(jd){
/* Reference:  The 1992 Astronomical Almanac, page B6. */
	//UT   = Frac(jd + 0.5);
	UT 		= Frac(jd + 0.5)
	jd   	= jd - UT
	TU   	= (jd - 2451545.0)/36525
	GMST 	= 24110.54841 + TU * (8640184.812866 + TU * (0.093104 - TU * 6.2E-6))
	//GMST 	= Modulus(GMST + secday*omega_E*UT,secday);
	GMST 	= (GMST + secday*omega_E*UT)%secday 

	return( twopi * GMST/secday );
} /*Function ThetaG_JD*/


function fmod2p(arg) {
	ret_val  	= arg
	i 			= parseInt(ret_val / twopi)
	ret_val   -= i * twopi;

	if (ret_val < 0) {
		ret_val += twopi;
	}

	return ret_val
}

/* Four-quadrant arctan function */
function AcTan(sinx, cosx) {
	if (cosx == 0) {
		if (sinx > 0) {
			return pio2
		} else {
			return x3pio2
		}
	} else {
		if (cosx > 0) {
			if (sinx > 0) {
				return Math.atan(sinx / cosx)
			} else {
				return twopi + Math.atan(sinx / cosx)
			}
		} else {
			return pi + Math.atan(sinx / cosx)
		}
	}
}

function magnitude(pos) {
	return Math.sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z)
} /*Procedure Magnitude*/


function geodetic2cartesian(lat,lng,alt, jd) {
/* Reference:  The 1992 Astronomical Almanac, page K11. */
	theta 		= fmod2p(thetag_JD(jd) + lng)/*LMST*/
	c 			= 1/Math.sqrt(1 + eflat*(eflat - 2)*Math.pow(Math.sin(lat),2))
	sq 			= Math.pow(1 - eflat, 2)*c
	achcp 		= (xkmper*c + alt)*Math.cos(lat)
	x 			= achcp*Math.cos(theta)/*kilometers*/
	y 			= achcp*Math.sin(theta)
	z 			= (xkmper*sq + alt)*Math.sin(lat)
	//magnitude(obs_pos)

	var result 	= {
		x 	: x,
		y 	: y,
		z 	: z,
		t 	: theta
	}
	return result	
}

/*
function cartesian2geodetic(satpos,ts)
{
	clConv 	= new Clock()
	jd 		= clConv.JulianDay(ts)
	// Reference:  The 1992 Astronomical Almanac, page K12.
	theta 	= Math.atan2(satpos.y,satpos.x)//radians
	lon 	= fmod2p(theta - thetag_JD(jd))//radians
	r 		= Math.sqrt(satpos.x*satpos.x + satpos.y*satpos.y)
	e2 		= eflat*(2 - eflat)
	lat 	= Math.atan2(satpos.z,r)//radians

	do
	{
		phi 	= lat
		c 		= 1/Math.sqrt(1 - e2*Math.sin(phi)*Math.sin(phi))
		lat 	= Math.atan2(satpos.z + xkmper*c*e2*Math.sin(phi),r)
	}
	while(Math.abs(lat - phi) >= 1e-10)

	alt = r/Math.cos(lat) - xkmper*c//kilometers

	if( lat > pio2 ) lat -= twopi

	var result 		=
	{
		lat		: lat,
		lng 	: lon,
		alt 	: alt,
	}
	return result
}*/ /*Procedure Calculate_LatLonAlt*/

//Конвертация декартовых координат в географические
function cartesian2geodetic(satpos,ts) {
	theta 	= AcTan(satpos.y,satpos.x)
	lng 	= fmod2p(theta - thetag_JD(JulianDay_calc(ts)))
	r 		= Math.sqrt(satpos.x*satpos.x + satpos.y*satpos.y)
	e2 		= eflat * (2 - eflat)
	lat 	= AcTan(satpos.z, r)

	/*console.log('tsince = '+app.prop.t)
	console.log('JD = '+JulianDay_calc(ts))
	console.log('theta = '+theta)
	console.log('r = '+r)
	console.log('e2 = '+e2)
	console.log('thetag = '+thetag_JD(JulianDay_calc(ts)))
	console.log('lon = '+lng)
	console.log('theta - thetag = '+(theta-thetag_JD(JulianDay_calc(ts))))*/

	do {
		phi 	= lat
		sinPhi 	= Math.sin(phi)
		c 		= 1 / (Math.sqrt(1 - e2 * (sinPhi*sinPhi)))
		lat 	= AcTan(satpos.z + xkmper * c * e2 * sinPhi, r)
	} while (Math.abs(lat - phi) >= 1e-10)

	alt = r / Math.cos(lat) - xkmper * c 

	if (lat > pio2) {
		lat -= twopi
	}

	satpos.lat = lat/torad
	satpos.lng = lng/torad
	satpos.alt = alt

	return satpos
	/* Reference:  The 1992 Astronomical Almanac, page K12.

	$geodetic->theta = Predict_Math::AcTan($pos->y, $pos->x); //radians
	$geodetic->lon = Predict_Math::FMod2p($geodetic->theta - Predict_Time::ThetaG_JD($_time)); //radian
	$r = sqrt(($pos->x * $pos->x) + ($pos->y * $pos->y));
	$e2 = Predict::__f * (2 - Predict::__f);
	$geodetic->lat = Predict_Math::AcTan($pos->z, $r); //radians

	do {
	$phi    = $geodetic->lat;
	$sinPhi = sin($phi);
	$c      = 1 / sqrt(1 - $e2 * ($sinPhi * $sinPhi));
	$geodetic->lat = Predict_Math::AcTan($pos->z + Predict::xkmper * $c * $e2 * $sinPhi, $r);
	} while (abs($geodetic->lat - $phi) >= 1E-10);

	$geodetic->alt = $r / cos($geodetic->lat) - Predict::xkmper * $c;//kilometers

	if ($geodetic->lat > Predict::pio2) {
	$geodetic->lat -= Predict::twopi;
	}
	 */

	/*
	var gmst 	= GMST_calc(ts)
	var lst		= gmst*15;

	var r 		= Math.sqrt(satpos.x*satpos.x+satpos.y*satpos.y)
	var lng 	= Math.atan2(satpos.y,satpos.x)/torad - lst
	if(lng>360)
	{
		lng 	= lng%360
	}
	if(lng<0)
	{
		lng 	= lng%360+360
	}    
	if(lng>180)
	{
		lng 	= lng-360
	}

	var lat 	= Math.atan2(satpos.z,r)
	var e2 		= eflat*(2-eflat)
	var tmp_lat = 0

	do
	{
		tmp_lat 	= lat
		var sin_lat = Math.sin(tmp_lat)
		var c 		= 1/Math.sqrt(1-e2*sin_lat*sin_lat)
		lat 		= Math.atan2(satpos.z+xkmper*c*e2*(Math.sin(tmp_lat)),r)
	}
	while(Math.abs(lat-tmp_lat)>0.0001)

	var alt 		= r/Math.cos(lat)-xkmper*c

	lat 			/= torad
	
	satpos.lat 		= lat
	satpos.lng 		= lng
	satpos.alt 		= alt

	return satpos
	*/
}

//Перевод координат модели SGP/SDP4 в киллометры и км/с
function cartesian2km(satpos) {
	xkm		= (satpos.x*xkmper)
	ykm 	= (satpos.y*xkmper)
	zkm 	= (satpos.z*xkmper)
	vxkm	= (satpos.xdot*xkmper/60)
	vykm 	= (satpos.ydot*xkmper/60)
	vzkm 	= (satpos.zdot*xkmper/60)

	var result 	= {
		x	: xkm,
		y 	: ykm,
		z 	: zkm,
		w 	: Math.sqrt(xkm*xkm+ykm*ykm+zkm*zkm),
		vx 	: vxkm,
		vy 	: vykm,
		vz 	: vzkm,
		vw 	: Math.sqrt(vxkm*vxkm+vykm*vykm+vzkm*vzkm)
	}
	return result
}

//Вычисление значений азимута и угла места
function observer_calc(eslat, eslong, satpos, ts) {
	eslat 			*= torad
	eslong 			*= torad
	jd 				= JulianDay_calc(ts)

	escartesian 	= geodetic2cartesian(eslat, eslong, 0, jd)

	diff = {
		x 	: satpos.x - escartesian.x,
		y 	: satpos.y - escartesian.y,
		z 	: satpos.z - escartesian.z		
	}

	geotheta 		= Math.atan2(escartesian.y,escartesian.x)/*radians*/


	sin_lat 		= Math.sin(eslat)
	cos_lat 		= Math.cos(eslat)
	sin_theta 		= Math.sin(geotheta)
	cos_theta 		= Math.cos(geotheta)

	top_s 	= sin_lat * cos_theta * diff.x + sin_lat * sin_theta * diff.y - cos_lat * diff.z
	top_e 	= -sin_theta * diff.x + cos_theta * diff.y
	top_z 	= cos_lat * cos_theta * diff.x + cos_lat * sin_theta * diff.y + sin_lat * diff.z

	posw  	= magnitude(diff)

	azim 	= Math.atan(-top_e/top_s) /*Azimuth*/
	if( top_s > 0 ) {
		azim = azim + pi
	}
	if( azim < 0 ) {
		azim = azim + twopi
	}
	el = Math.asin(top_z/posw)

	var result 	= {
		azimuth		: azim/torad,
		elevation 	: el/torad,
	}
	return result
} /*Procedure Calculate_Obs*/

function controlURL() {
	console.log('http://sandbox.roundeasy.ru/predictPHP/?tsince='+app.prop.t+'&JD='+JulianDay_calc(curTimeStampUTC()))
}