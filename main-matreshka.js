var Application = Class({
	'extends': Matreshka,
	constructor: function() {
		this.mapSelector = 'map-canvas'
		this.satsSelector = 'sats'

		this.reset()

		this.startTime	 	= curTimeStampUTC()-86400*1000
		this.stopTime 		= curTimeStampUTC()+86400*1000

		//Связи с элементами управления
		this.bindNode( 'tleFile', '#file' )
		this.bindNode( 'satName', '#sats' )
		//Кнопки
		this.bindNode( 'chooseMarkerBtn', '#chooseMarker' )
		this.bindNode( 'choosePosWindowBtn', '#choosePosWindowBtn' )
		this.bindNode( 'controlTitleBtn', '#controlTitle' )
		this.bindNode( 'showSatBtn', '#showSat' )
		this.bindNode( 'showMarkerBtn', '#showMarker' )
		
		//Вывод
		this.bindNode( 'obsLat', '#obsLat', {
			setValue: function(val) {
				this.innerHTML = val
			}
		})
		this.bindNode( 'obsLng', '#obsLng', {
			setValue: function(val) {
				this.innerHTML = val
			}
		})

		//Установка точки наблюдения
		this.on( 'click::choosePosWindowBtn', function() {
			if (this.choosePosMode) {
				lat = parseFloat(document.getElementById("window-obs-lat").value)
				lng = parseFloat(document.getElementById("window-obs-lng").value)
				this.setObserver(lat, lng)
			}
		})
		//Скрытие раскрытие панели управления
		this.on( 'click::controlTitleBtn', function() {

			windowMode   		= document.getElementById("controls").style.display
			windowMode 			= windowMode=='none'?'block':'none'
			document.getElementById("controls").style.display = windowMode			
		})
		//Отображение спутника с точки на Земле
		this.on( 'click::showSatBtn', function() {
			this.showSatFromEarth()
		})
		//Отображение точки на Земле со спутника
		this.on( 'click::showMarkerBtn', function() {
			this.showEarthFromSat()
		})
		//Кнопка выбора точка наблоюдения
		this.on( 'click::chooseMarkerBtn', function ()	{
  
			this.choosePosMode = !this.choosePosMode
			windowMode   		= this.choosePosMode?'block':'none'
			document.getElementById("choosePosWindow").style.display = windowMode
		})
		//Выбор TLE-файла
		this.on( 'change:tleFile', function() {
			if (this.tleFile == 0) return

			this.earth.clock.shouldAnimate 	= false
			this.tick 						= 0

			this.readTLE(this.tleFile)
		})
		//Выбор спутника
		this.on( 'change:satName', function() {
			if (this.satName == 0) return

			this.earth.clock.shouldAnimate 	= false
			this.tick 						= 0
			document.getElementById("obs").style.display = 'none' 
			this.preCalcSat(this.satName)
		})
	},
	getColor: function() {
		if (typeof this.lastColor === 'undefined') {
			this.lastColor = 0
		}

		switch (this.lastColor) {
			case 0:
				color = Cesium.Color.GREEN.withAlpha(0.2)
				break
			case 1:
				color = Cesium.Color.BLUE.withAlpha(0.2)
				break
			case 2:
				color = Cesium.Color.YELLOW.withAlpha(0.2)
				break
			case 3:
				color = Cesium.Color.MEDIUMVIOLETRED.withAlpha(0.2)
				break
			default:
				color = Cesium.Color.GREEN.withAlpha(0.2)
				this.lastColor = 0
		}
		this.lastColor++

		return color
	},
	preCalcSat: function(satName) {
		if (typeof this.selSats[satName] === 'undefined')
		{
			//Получаем строки TLE файла для текущего спутника
			tle_line1 				= this.sats[satName].line1
			tle_line2 				= this.sats[satName].line2
			//Выделяем орбитальные параметры из TLE в переменнную orbital_elements
			orbital_elements 		= TLEParse(tle_line1, tle_line2)
			//Инициализация алгоритма SGP4/SDP4
			this.selSats[satName] 	= new SGDP4(orbital_elements)		
		}
		this.satName = satName
		this.calc24(satName)
		this.calcAll()
	},
	setObserver: function(lat, lng) {
		if (this.earthInitFlag) {
			this.earth.entities.remove(this.earthEntity)

			var pinBuilder = new Cesium.PinBuilder();
			this.earthEntity = this.earth.entities.add({
				name : 'Точка наблюдения',
				position : Cesium.Cartesian3.fromDegrees(lng, lat),
				billboard : {
					image : pinBuilder.fromText('*', Cesium.Color.BLACK, 48).toDataURL(),
					verticalOrigin : Cesium.VerticalOrigin.BOTTOM
				}
			})
		}

		this.observer 		= {}
		this.observer.lng 	= round(lng,2)
		this.observer.lat 	= round(lat,2)

		this.obsLat 		= this.observer.lat
		this.obsLng 		= this.observer.lng
		this.choosePosMode 	= false	

		document.getElementById("choosePosWindow").style.display = 'none' 
		document.getElementById("obs").style.display 	= 'block' 

		this.calcAll()
		//this.calcObserver()
	},
	calcObserver: function(satName, ts) {
		if (typeof this.observer.lng !== 'undefined' &&
			typeof this.observer.lat !== 'undefined') {
			var aztRes 			= observer_calc(this.observer.lat, 
										this.observer.lng, this.satPos[satName], ts)

			this.satPos[satName].azimuth 	= round(aztRes.azimuth,2)
			this.satPos[satName].elevation 	= round(aztRes.elevation,2)
		}
	},
	initEarth: function() {
		this.earth 			= new Cesium.Viewer('cesiumContainer', {
		//Use OpenStreetMaps
			imageryProvider : new Cesium.OpenStreetMapImageryProvider({
				url : "lib/Cesium/Assets/Textures/NaturalEarthII/tilemapresource.xml"
			}),
		})


		this.earthInitFlag 	= true

		this.earth.canvas.setAttribute('tabindex', '0')
		this.earth.scene.globe.enableLighting 	= true
		this.earth.clock.shouldAnimate 			= false
		this.tick 								= 0

		this.flags = {
			looking : false,
			moveForward : false,
			moveBackward : false,
			moveUp : false,
			moveDown : false,
			moveLeft : false,
			moveRight : false,
			rotateLeft : false,
			rotateRight : false
		}
		//this.earth.scene.globe.enableLighting = true
	},
	readTLE: function(file)	{
		var sats 	= {}
		parent 		= this
		this.sats 	= sats
		var rawFile = new XMLHttpRequest();
		var sel 	= document.getElementById(this.satsSelector);

		rawFile.open("GET", file, true);
		rawFile.onreadystatechange = function () {
			if(rawFile.readyState === 4) {
				if(rawFile.status === 200 || rawFile.status == 0) {
					var allText 	= rawFile.responseText
					var lines 		= allText.split('\n')
					var curLine 	= '',
						sat_name	= '',
						line1 		= '',
						line2 		= '',
						sat 		= new Array()

					sel.innerHTML = '<option value="0" selected>Выберите спутник</option>'
					for(var line = 0, num = 1; line < lines.length; line++,num++)
					{
						var curLine = lines[line]
						if (num == 1) {
							sat_name = curLine
						} else if (num == 2) {
							line1 = curLine
						} else if (num == 3) {
							sat_name 	= sat_name.trim()
							line2 		= curLine
							sat = {
								title 	: sat_name,
								line1 	: line1,
								line2 	: line2
							}
							sats[sat_name] = sat
							var option 	= document.createElement("option")
							option.text = sat_name
							sel.add(option)
							num = 0
						}
					}
					parent.sats = sats
				}
			}
		}
		rawFile.send(null)
	},
	delSat: function(index) {
		//Удаляем из списка
		delete this.selSats[index]
		delete this.satPos[index]
		//Удаляем с модели
		this.earth.entities.remove(this.satToEarthEntities[index])
		this.earth.entities.remove(this.bodyEntities[index])
		this.earth.entities.remove(this.coverEntities[index])
		delete(this.bodyEntities[index])
		//Пересчет
		this.calcAll()
	},
	selSat: function(index) {
		this.satName = index
	},
	calcAll: function() {
		if (Object.keys(app.selSats).length > 0) {
			document.getElementById("sat").style.display 	= 'block' 
		} else {
			document.getElementById("sat").style.display 	= 'none' 
		}

		document.getElementById("sat-table").innerHTML = ''
		for(var index in this.selSats) { 
			if (this.selSats.hasOwnProperty(index)) {
				this.calc(index)

				if (typeof this.satPos[index].azimuth == 'undefined') this.satPos[index].azimuth = '-'
				if (typeof this.satPos[index].elevation == 'undefined') this.satPos[index].elevation = '-'

				var tr = document.createElement("tr")
				tr.innerHTML = '<td>'+index+'</td>'+
					'<td>'+this.satPos[index].lat+'</td>'+
					'<td>'+this.satPos[index].lng+'</td>'+
					'<td>'+this.satPos[index].alt+'</td>'+
					'<td>'+this.satPos[index].azimuth+'</td>'+
					'<td>'+this.satPos[index].elevation+'</td>'+
					'<td><a href="#" onclick="app.delSat(\''+index+'\'); return false;">Удалить</a>'+' | '+
					'<a href="#" onclick="app.selSat(\''+index+'\'); return false;">Выбрать</a></td>'
				document.getElementById("sat-table").appendChild(tr)
			}
		}
	},
	calc: function(satName, setTimeStamp)	{
		//Получаем текущее время или устанавливаем из переменной
		if (typeof setTimeStamp !== 'undefined') {
			var curTimeStamp 	= setTimeStamp
		} else {
			var curTimeStamp 	= Cesium.JulianDate.toDate(this.earth.clock.currentTime).getTime()
		}
		//Расчет текущей позиции спутника
		var propPos 			= this.selSats[satName].calc(curTimeStamp)
		//Перевод единиц радиуса Земли в км (скорости и позиции)
		this.satPos[satName]	= cartesian2km(propPos)
		//Перевод в географические координаты из декартовой системы
		cartesian2geodetic(this.satPos[satName], curTimeStamp)
		//Вычисление зоны покрытия
		this.satPos[satName] 	= footprint_calc(this.satPos[satName])

		//Округление полученных значений географических координат
		this.satPos[satName].lat 	= round(this.satPos[satName].lat,2)
		this.satPos[satName].lng 	= round(this.satPos[satName].lng,2)
		this.satPos[satName].alt 	= round(this.satPos[satName].alt,0)

		//Вычисление азимута и угла места для заданного спутника и времени
		this.calcObserver(satName, curTimeStamp)

		return this.satPos[satName]
	},
	reset:  function()	{
		this.tleFile 			= 0
		this.satName 			= 0

		this.selSats 			= {}
		this.satPos 			= {}

		this.bodyEntities 		= {}
		this.coverEntities 		= {}
		this.satToEarthEntities = {}

		this.choosePosMode 		= false
		this.earthInitFlag 		= false
		this.followSat 			= false

		this.observer 			= {}

		document.getElementById("obs").style.display 	= 'none' 
		document.getElementById("sat").style.display 	= 'none' 

		this.initEarth()
	},
	calc_orbit: function(satName)	{
		var curTimeStamp 	= curTimeStampUTC()
		var property 		= new Cesium.SampledPositionProperty()
		var surfaceProperty = new Cesium.SampledPositionProperty()
		if (this.earthInitFlag) {
			for (ts = curTimeStamp-86400*1000; ts <= curTimeStamp+86400*1000; ts += 3600*1000) {
				var time 			= Cesium.JulianDate.fromDate(dateFromTimestamp(ts))
				var coors 			= this.calc(satName, ts)
				var position 		= Cesium.Cartesian3.fromDegrees(coors.lng, coors.lat, coors.alt*1000)
				var surfacePosition = Cesium.Cartesian3.fromDegrees(coors.lng, coors.lat, 0)

				if (ts == curTimeStamp) {
					var coverage 		= coverage_calc(coors, '3D')
				}				
				property.addSample(time, position)
				surfaceProperty.addSample(time, surfacePosition)
			}
			return {sat: property, surface: surfaceProperty, coverage: coverage}
			//polyline.positions 		= new Cesium.ConstantProperty(coordinatesArr)
		}
	},
	calc24: function(satName) {
		start 							= Cesium.JulianDate.fromDate(dateFromTimestamp(this.startTime))
		stop 							= Cesium.JulianDate.fromDate(dateFromTimestamp(this.stopTime))
		current 						= Cesium.JulianDate.fromDate(dateFromTimestamp(curTimeStampUTC()))
		this.earth.clock.startTime 		= start.clone()
		this.earth.clock.stopTime 		= stop.clone()
		this.earth.clock.currentTime 	= current.clone()
		this.earth.clock.clockRange 	= Cesium.ClockRange.LOOP_STOP
		this.earth.timeline.zoomTo(start, stop)

		//Цвет
		if (typeof this.selSats[satName].useColor === 'undefined') {
			this.selSats[satName].useColor = this.getColor()
		}
		//this.earth.entities.remove(this.bodyEntities[satName])
		//Compute the entity position property.
		var positionArr = this.calc_orbit(satName)

		//Положение модели спутника
		var heading = Cesium.Math.toRadians(0)
		var pitch = Cesium.Math.toRadians(0)
		var roll = Cesium.Math.toRadians(90)
		var orientation = Cesium.Transforms.headingPitchRollQuaternion(positionArr.sat, heading, pitch, roll)
				
		positionArr.sat.setInterpolationOptions({
			interpolationDegree : 5,
			interpolationAlgorithm : Cesium.LagrangePolynomialApproximation
		})
		positionArr.surface.setInterpolationOptions({
			interpolationDegree : 5,
			interpolationAlgorithm : Cesium.LagrangePolynomialApproximation
		})
		//Actually create the entity
		this.bodyEntities[satName] = this.earth.entities.add({
			name : satName,
			//Use our computed positions
			position : positionArr.sat,
			//Automatically compute orientation based on position movement.
			orientation : new Cesium.VelocityOrientationProperty(positionArr.sat),
			//Set the entity availability to the same interval as the simulation time.
			availability : new Cesium.TimeIntervalCollection([new Cesium.TimeInterval({
				start : start,
				stop : stop
			})]),
			//Load the Cesium plane model to represent the entity
			model : {
				uri : 'lib/sat.gltf',
				minimumPixelSize : 64
			},
		})

		this.coverEntities[satName] = this.earth.entities.add({
			position: positionArr.surface,
			name : satName,
			//Set the entity availability to the same interval as the simulation time.
			availability : new Cesium.TimeIntervalCollection([new Cesium.TimeInterval({
				start : start,
				stop : stop
			})]),
			ellipse : {
				semiMinorAxis : positionArr.coverage.semiMinor/2,
				semiMajorAxis : positionArr.coverage.semiMajor/2,
				material : new Cesium.GridMaterialProperty({
					color : this.selSats[satName].useColor,
				}),//this.selSats[satName].useColor,
				outline : true,
				outlineColor : Cesium.Color.BLACK,
			}
		})	
	},
	showSatFromEarth: function() {
		if (typeof this.satPos[this.satName] == 'undefined' ||
			typeof this.satPos[this.satName].lng == 'undefined' || 
			typeof this.satPos[this.satName].lat == 'undefined') {	
			alert('Выберите спутник')
			return false
		}
		if (typeof this.observer.lng == 'undefined' ||
			typeof this.observer.lat == 'undefined') {
			earthPos 	= Cesium.Cartesian3.fromDegrees(this.satPos[this.satName].lng, this.satPos[this.satName].lat, 0)
		} else {
			earthPos 	= Cesium.Cartesian3.fromDegrees(this.observer.lng, this.observer.lat, 0)
		}

		this.earth.camera.setView({
			position : earthPos,
			heading : Cesium.Math.toRadians(0.0),
			pitch : Cesium.Math.toRadians(0),
			roll : 0.0
		})
		satPosition	= Cesium.Cartesian3.fromDegrees(this.satPos[this.satName].lng, this.satPos[this.satName].lat, this.satPos[this.satName].alt*1000)
		this.earth.camera.lookAt(satPosition)		
	},
	showEarthFromSat: function() {
		if (typeof this.satPos[this.satName] == 'undefined' ||
			typeof this.satPos[this.satName].lng == 'undefined' || 
			typeof this.satPos[this.satName].lat == 'undefined') {	
			alert('Выберите спутник')
			return false
		}
		if (typeof this.observer.lng == 'undefined' ||
			typeof this.observer.lat == 'undefined') {
			earthPos 	= Cesium.Cartesian3.fromDegrees(this.satPos[this.satName].lng, this.satPos[this.satName].lat, 0)
		} else {
			earthPos 	= Cesium.Cartesian3.fromDegrees(this.observer.lng, this.observer.lat, 0)
		}

		satPosition	= Cesium.Cartesian3.fromDegrees(this.satPos[this.satName].lng, this.satPos[this.satName].lat, this.satPos[this.satName].alt*1000)
		this.earth.camera.setView({
			position : satPosition,
			heading : Cesium.Math.toRadians(90),
			pitch : Cesium.Math.toRadians(0),
			roll : 0.0
		})
		this.earth.camera.lookAt(earthPos)		
	}
})

var app = new Application()

app.earth.clock.onTick.addEventListener(function(clock) {
	if (app.earth.clock.shouldAnimate) {
		if (app.tick>=10)
		{
			app.tick = 0
			app.calcAll()			
		} else {
			app.tick++
		}
	}
	// Change movement speed based on the distance of the camera to the surface of the ellipsoid.
	var cameraHeight = app.earth.scene.globe.ellipsoid.cartesianToCartographic(app.earth.camera.position).height
	var moveRate = cameraHeight / 100.0

	if (app.flags.moveForward) {
		app.earth.camera.moveForward(moveRate)
	}
	if (app.flags.moveBackward) {
		app.earth.camera.moveBackward(moveRate)
	}
	if (app.flags.moveUp) {
		app.earth.camera.moveUp(moveRate)
	}
	if (app.flags.moveDown) {
		app.earth.camera.moveDown(moveRate)
	}
	if (app.flags.moveLeft) {
		app.earth.camera.moveLeft(moveRate)
	}
	if (app.flags.moveRight) {
		app.earth.camera.moveRight(moveRate)
	}
	if (app.flags.rotateLeft) {
		app.earth.camera.rotateLeft()
	}
	if (app.flags.rotateRight) {
		app.earth.camera.rotateRight()
	}
});	
document.addEventListener('keydown', function(e) {
	var flagName = getFlagForKeyCode(e.keyCode)
	if (typeof flagName !== 'undefined') {
		app.flags[flagName] = true
	}
}, false)
document.addEventListener('keyup', function(e) {
	var flagName = getFlagForKeyCode(e.keyCode)
	if (typeof flagName !== 'undefined') {
		app.flags[flagName] = false
	}
}, false)