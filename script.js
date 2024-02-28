var lsib = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017"),
    terra = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE"),
    visVPD = {"opacity":1,"bands":["vpd_sc"],"min":0,"max":30,"palette":["e5f8d2","7ddc1f","ffff00","ffa500","ff2500"]},
    visHDW = {"opacity":1,"bands":["HDW"],"min":0,"max":8,"palette":["e5f8d2","7ddc1f","ffff00","ffa500","ff2500"]},
    visvs = {"opacity":1,"bands":["vs"],"min":0,"max":350,"palette":["e5f8d2","7ddc1f","ffff00","ffa500","ff2500"]},
    visrh = {"opacity":1,"bands":["rh"],"min":0,"max":120,"palette":["ff2500","ffa500","ffff00","7ddc1f","e5f8d2"]},
    vispr = {"opacity":1,"bands":["sc_pr"],"min":0,"max":80,"palette":["e5f8d2","7ddc1f","ffff00","ffa500","ff2500"]},
    vislfdi = {"opacity":1,"bands":["LFDI"],"min":20,"max":70,"palette":["e5f8d2","7ddc1f","ffff00","ffa500","ff2500"]};

//---------------------------------------------------------General Functions

//function to remove layers from map
function removelay(){
  var lay = Map.layers().get(0);
  if(lay){
  Map.remove(lay)}

  else{print('layer missing')}
}

//Vapor pressure deficit scaling
function scaleVPD(img){
  var scVPD = img.expression(
    'vpd*0.01', 
    {
     'vpd': img.select('vpd')
    }
    ).rename('vpd_sc');
  
  return img.addBands(scVPD);
  } 
  
//HDW index calculation (wirh wind scaling)
function addHDW(img){
  var HDW = img.expression(
    'vpd*u*0.01', 
    {
     'vpd': img.select('vpd_sc'),
     'u': img.select('vs')
    }
    ).rename('HDW');
  
  return img.addBands(HDW);
}

//Mean temperature estimation (with temperature scaling)
function addT(img){
  var T = img.expression(
    '((Tmax*0.1)+(Tmin*0.1))/2', 
    {
     'Tmin': img.select('tmmn'),
     'Tmax': img.select('tmmx')
    }
    ).rename('Temp');
  
  return img.addBands(T);
  }

//Calculation of saturated vapor pressure from Tetens equation
function addes(img){
  var es = img.expression(
    '0.61078*exp((17.269 *T)/(237.3+T))', 
    {
     'T': img.select('Temp')
    }
    ).rename('es');
  
  return img.addBands(es);
  }

//Calculation of relative humidity from VPD and saturated vapor pressure
function addrh(img){
  var rh = img.expression(
    '(1-(vpd/es))*100', 
    {
     'vpd': img.select('vpd_sc'),
     'es': img.select('es')
    }
    ).rename('rh');
  
  return img.addBands(rh);
  }

//Calculation of BI parameter of LFDI index 
function addbi(img){
  var BI = img.expression(
    '(T-35)-((35-T)/30)+((100-RH)*0.37)+30', 
    {
     'T': img.select('Temp'),
     'RH': img.select('rh')
    }
    ).rename('bi');
  
  return img.addBands(BI);
  }

//Wind speed conversion to km/h and scaling 
function scaleWind(img){
  var ws = img.expression(
    'v*0.01*3.6', 
    {
     'v': img.select('vs')
    }
    ).rename('ws');
  
  return img.addBands(ws);
  } 
  
//Calculation of WF parameter of LFDI index 
function addwf(img){
  var WF = img.expression(
    '-0.0000227*(ws**4)+0.0026348*(ws**3)-0.09087*(ws**2)+1.65*ws+0.2', 
    {
     'ws': img.select('ws')
    }
    ).rename('wf');
  
  return img.addBands(WF);
  }

//Scale extreme rain accumulation values  
function scalepr(img){
 var unsc_pr = img.select('pr').rename('sc_pr');
 var sc_pr = unsc_pr
  .where(unsc_pr.gt(100),100);
  
  return img.addBands(sc_pr);
}

//Calculation of RDF parameter of LFDI index 
function addrcf(img){
  var rcf = img.expression(
    '0.62-(0.0342*p)+0.000609*(p**2)-0.000004*(p**3)+0.1761*10-0.01141*(10**2)+0.000279*(10**3)', 
    {
     'p': img.select('sc_pr')
    }
    ).rename('rcf');
  
  return img.addBands(rcf);
  }
  
  //Calculation of LFDI index   
function addLFDI(img){
  var LFDI = img.expression(
    '((bi+wf)*rcf)', 
    {
     'wf': img.select('wf'),
     'bi': img.select('bi'),
     'rcf': img.select('rcf'),
    }
    ).rename('LFDI');
  
  return img.addBands(LFDI);
  }  

//Clip to area of interest
function clipAOI(img){
  return img.clip(Greece);
}

//Area of interest definition
var Greece = lsib.filter(ee.Filter.eq('country_na','Greece'));

//Dataset subseting and caluculation of parameters
var terra_subset = terra.filter(ee.Filter.calendarRange(2018,2023,'year'))
  .filter(ee.Filter.calendarRange(5,9,'month'))
  .map(scaleVPD)
  .map(addHDW)
  .map(addT)
  .map(addes)
  .map(addrh)
  .map(addbi)
  .map(scalepr)
  .map(scaleWind)
  .map(addwf)
  .map(addrcf)
  .map(addLFDI)
  .map(clipAOI);
  

//---------------------------------Results Functions (all follow the same patern)
function may(){
  removelay();
  removelay();
  removelay();
  
  //Calculation of mean wildfire risk indices
  var may_mean = terra_subset.filter(ee.Filter.calendarRange(5,5,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.mean());
  
  //Calculation of correlation of LFDI and HDW
  var may_corr = terra_subset.filter(ee.Filter.calendarRange(5,5,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.pearsonsCorrelation());
  
  //Add results to map
  Map.addLayer(may_corr,{},'Βαθμός συσχέτισης LFDI-HDW');
  Map.addLayer(may_mean.select('HDW_mean').rename('HDW'),visHDW,'Μέσος HDW');
  Map.addLayer(may_mean.select('LFDI_mean').rename('LFDI'),vislfdi,'Μέσος LFDI');
  
  // Export.image.toDrive({
  //   image: may_mean.select('LFDI_mean'),
  //   description: 'LFDI_mean_may',
  //   scale: 4500,
  //   region: Greece,
  //   maxPixels: 1e13,
  // });
}

function june(){
  removelay();
  removelay();
  removelay();
  var may_mean = terra_subset.filter(ee.Filter.calendarRange(6,6,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.mean());
  var may_corr = terra_subset.filter(ee.Filter.calendarRange(6,6,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.pearsonsCorrelation());

  Map.addLayer(may_corr,{},'Βαθμός συσχέτισης LFDI-HDW');
  Map.addLayer(may_mean.select('HDW_mean').rename('HDW'),visHDW,'Μέσος HDW');
  Map.addLayer(may_mean.select('LFDI_mean').rename('LFDI'),vislfdi,'Μέσος LFDI');
  
  // Export.image.toDrive({
  //   image: may_mean.select('LFDI_mean'),
  //   description: 'LFDI_mean_june',
  //   scale: 4500,
  //   region: Greece,
  //   maxPixels: 1e13,
  // });
}

function july(){
  removelay();
  removelay();
  removelay();
  var may_mean = terra_subset.filter(ee.Filter.calendarRange(7,7,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.mean());
  var may_corr = terra_subset.filter(ee.Filter.calendarRange(7,7,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.pearsonsCorrelation());

  Map.addLayer(may_corr,{},'Βαθμός συσχέτισης LFDI-HDW');
  Map.addLayer(may_mean.select('HDW_mean').rename('HDW'),visHDW,'Μέσος HDW');
  Map.addLayer(may_mean.select('LFDI_mean').rename('LFDI'),vislfdi,'Μέσος LFDI');
  
  // Export.image.toDrive({
  //   image: may_mean.select('LFDI_mean'),
  //   description: 'LFDI_mean_july',
  //   scale: 4500,
  //   region: Greece,
  //   maxPixels: 1e13,
  // });
}

function august(){
  removelay();
  removelay();
  removelay();
  var may_mean = terra_subset.filter(ee.Filter.calendarRange(8,8,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.mean());
  var may_corr = terra_subset.filter(ee.Filter.calendarRange(8,8,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.pearsonsCorrelation());

  Map.addLayer(may_corr,{},'Βαθμός συσχέτισης LFDI-HDW');
  Map.addLayer(may_mean.select('HDW_mean').rename('HDW'),visHDW,'Μέσος HDW');
  Map.addLayer(may_mean.select('LFDI_mean').rename('LFDI'),vislfdi,'Μέσος LFDI');

  // Export.image.toDrive({
  //   image: may_mean.select('LFDI_mean'),
  //   description: 'LFDI_mean_august',
  //   scale: 4500,
  //   region: Greece,
  //   maxPixels: 1e13,
  // });
  
}

function september(){
  removelay();
  removelay();
  removelay();
  var may_mean = terra_subset.filter(ee.Filter.calendarRange(9,9,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.mean());
  var may_corr = terra_subset.filter(ee.Filter.calendarRange(9,9,'month')).select(['LFDI','HDW']).reduce(ee.Reducer.pearsonsCorrelation());

  Map.addLayer(may_corr,{},'Βαθμός συσχέτισης LFDI-HDW');
  Map.addLayer(may_mean.select('HDW_mean').rename('HDW'),visHDW,'Μέσος HDW');
  Map.addLayer(may_mean.select('LFDI_mean').rename('LFDI'),vislfdi,'Μέσος LFDI');
  
  // Export.image.toDrive({
  //   image: may_mean.select('LFDI_mean'),
  //   description: 'LFDI_mean_september',
  //   scale: 4500,
  //   region: Greece,
  //   maxPixels: 1e13,
  // });
}

//-------------------------UI Parameters
//Create a control panel
var ControlPanel = ui.Panel({
  style: {
    backgroundColor: 'white',
    border: '1px solid black',
    padding: '5px',
    width: '250px',
    position:'top-right'
  }
});

Map.add(ControlPanel);

//Activation buttons
var may_run = ui.Button({
  label: "Μαίος",
  onClick: may, 
  style: {
    stretch: "horizontal",
    height:'50px',
    fontWeight:'50px',
    Color:'#C11B17'
  }
});

var june_run = ui.Button({
  label: "Ιούνιος",
  onClick: june, 
  style: {
    stretch: "horizontal",
    height:'50px',
    fontWeight:'50px',
    Color:'#C11B17'
  }
});

var july_run = ui.Button({
  label: "Ιούλιος",
  onClick: july, 
  style: {
    stretch: "horizontal",
    height:'50px',
    fontWeight:'50px',
    Color:'#C11B17'
  }
});

var august_run = ui.Button({
  label: "Αύγουστος",
  onClick: august, 
  style: {
    stretch: "horizontal",
    height:'50px',
    fontWeight:'50px',
    Color:'#C11B17'
  }
});

var september_run = ui.Button({
  label: "Σεπτέμβριος",
  onClick: september, 
  style: {
    stretch: "horizontal",
    height:'50px',
    fontWeight:'50px',
    Color:'#C11B17'
  }
});

ControlPanel.add(may_run);
ControlPanel.add(june_run);
ControlPanel.add(july_run);
ControlPanel.add(august_run);
ControlPanel.add(september_run);
