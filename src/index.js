'use strict';

const { readFileSync, writeFileSync } = require('fs');
const { join } = require('path');

const { convertFileList } = require('brukerconverter');
const { fileListFromPath } = require('filelist-from');
const { gsd, optimizePeaks } = require('ml-gsd');
const { xyExtract } = require('ml-spectra-processing');
const {
  solventSuppression,
  xyAutoPeaksPicking,
  xyzAutoZonesPicking,
  xyzJResAnalyzer,
} = require('nmr-processing');
const { SpectrumGenerator } = require('spectrum-generator');

const path = 'C:\\Users\\alejo\\Documents\\BIOGUNE\\20S01360';
const pathToWrite = './';

const csvData = readFileSync('./src/annotationDB.csv', 'utf-8');
const database = getJSON(csvData);
const ROI = getROIs(database, 'urine', [
  { name: 'delta [ppm]', saveAs: 'delta' },
  { name: 'from [ppm]', saveAs: 'from' },
  { name: 'to [ppm]', saveAs: 'to' },
]);

let process2DOptions = {
  zonesPicking: {
    tolerances: [5, 100],
    nuclei: ['1H', '1H'],
    realTopDetection: false,
  },
  jResAnalyzer: { reference: 0, getZones: true },
};

let gsdOptions = {
  minMaxRatio: 0.01,
  broadRatio: 0.00025,
  smoothY: true,
  realTopDetection: true,
};

let converterOptions = {
  converter: { xy: true },
  filter: {
    experimentNumber: [11, 12],
    processingNumber: [1],
    ignoreFID: true,
    ignore2D: false,
  },
};

let optimizationOptions = {
  factorWidth: 1,
  factorLimits: 2,
  shape: {
    kind: 'gaussian',
  },
  optimization: {
    kind: 'lm',
    options: {
      maxInterations: 500,
      gradientDifference: 1e-3,
    },
  },
};

let alignmentOptions = {
  // reference peaks is the pattern to use only relative intensity import
  referencePeaks: [
    { x: 5.226, y: 1 },
    { x: 5.22, y: 1 },
  ],
  // the expected delta of reference signal,
  delta: 5.22,
  // the region to make the PP and search the reference signal
  fromTo: { from: 5.1, to: 5.4 },
};

async function main() {
  const fileList = fileListFromPath(path);
  const pdata = await convertFileList(fileList, converterOptions);

  const groups = {};
  for (let data of pdata) {
    if (!groups[data.name]) groups[data.name] = [];
    groups[data.name].push(data);
  }

  for (const group in groups) {
    let result = [];
    for (const data of groups[group]) {
      if (data.twoD) {
        result.push(process2D(data, process2DOptions));
      } else {
        result.push(process1D(data));
      }
    }
    writeFileSync(
      join(pathToWrite, `result_${result[0].name}.json`),
      JSON.stringify(result),
    );
  }
}

main();

function getROIs(data, target, mapping) {
  let rois = [];
  for (let roi of data) {
    if (roi[target]) {
      const targetData = roi[target];
      let tmp = {};
      if (Number.isFinite(targetData[mapping[0].name])) {
        for (let map of mapping) {
          tmp[map.saveAs] = targetData[map.name];
        }
        rois.push(tmp);
      }
    }
  }
  return rois;
}
function getJSON(data) {
  const lines = data.split('\n');
  let headers = [];
  let secondHeaders = [];
  let secondHeadersData = lines[1].split(',');
  const firstHeadersData = lines[0].split(',');
  for (let i = 0; i < firstHeadersData.length; i++) {
    const cell = firstHeadersData[i];
    if (cell.length > 0) {
      headers.push({
        fromIndex: i,
        name: cell,
      });
    }
  }

  for (let i = 1; i < headers.length - 1; i++) {
    headers[i - 1].toIndex = headers[i].fromIndex - 1;
  }

  for (let cell of secondHeadersData) {
    secondHeaders.push({
      name: cell.replace('\r', ''),
    });
  }

  const getHeaderName = (headers, index) => {
    for (const header of headers) {
      const { fromIndex, toIndex } = header;
      if (index >= fromIndex && index <= toIndex) return header.name;
    }
    return 'null';
  };

  const getValue = (str) => {
    if (str.length === 0) return str;
    return !isNaN(Number(str)) ? Number(str) : str;
  };

  let metabolites = [];
  for (let i = 2; i < lines.length; i++) {
    const cells = lines[i].split(',');
    let metabolite = {};
    for (let j = 0; j < secondHeaders.length; j++) {
      let headerName = getHeaderName(headers, j);
      if (!metabolite[headerName]) metabolite[headerName] = {};
      metabolite[headerName][secondHeaders[j].name] = getValue(
        cells[j].replace('\r', ''),
      );
    }
    metabolite = { ...metabolite.null, ...metabolite };
    delete metabolite.null;
    metabolites.push(metabolite);
  }

  return metabolites;
}

function align(input) {
  const { spectrum, referencePeaks, delta, fromTo } = input;

  const xyData = { x: spectrum.x, y: spectrum.re };
  const peaks = xyAutoPeaksPicking(xyData, {
    ...fromTo,
    optimize: true,
    shape: { kind: 'lorentzian' },
    groupingFactor: 2.5,
  });

  const marketPeaks = solventSuppression(
    peaks,
    [
      {
        delta,
        peaks: referencePeaks,
      },
    ],
    { markSolventPeaks: true },
  );

  if (peaks.length > 0) {
    const glucosePeaks = marketPeaks.filter((peak) => peak.kind === 'solvent');
    if (glucosePeaks.length < 1) {
      throw new Error('glucose peaks had not been found');
    }
    const shift =
      delta - glucosePeaks.reduce((a, b) => a + b.x, 0) / glucosePeaks.length;
    xyData.x.forEach((e, i, arr) => (arr[i] += shift));
  }

  return xyData;
}

function process1D(data) {
  let spectrum = data.spectra[0].data;
  if (spectrum.x[0] > spectrum.x[1]) {
    spectrum.x = spectrum.x.reverse();
    spectrum.re = spectrum.re.reverse();
  }

  const xyData = align({
    spectrum,
    ...alignmentOptions,
  });

  let peakOptimized = {
    name: data.source.name,
    expno: data.source.expno,
    fit: [],
  };
  for (let roi of ROI) {
    let experimental = xyExtract(xyData, { zones: [{ ...roi }] });

    let peaks = gsd(experimental, gsdOptions);
    let optPeaks = optimizePeaks(experimental, peaks, optimizationOptions);
    let spectrumGenerator = new SpectrumGenerator({ ...roi, nbPoints: 256 });

    for (let peak of optPeaks) {
      spectrumGenerator.addPeak({ ...peak });
    }

    peakOptimized.fit.push({
      roi,
      peaks,
      optimizedPeaks: optPeaks,
      experimental,
      fitted: spectrumGenerator.getSpectrum(),
    });
  }
  return peakOptimized;
}

function process2D(data, options) {
  const { zonesPicking, jResAnalyzer } = options;

  const frequency = data.meta.observeFrequency;
  const minMax = data.minMax;
  const zones = xyzAutoZonesPicking(minMax, {
    observedFrequencies: [frequency, frequency],
    ...zonesPicking,
  });

  const newZones = [];
  for (let zone of zones) {
    newZones.push(
      xyzJResAnalyzer(zone.signals, {
        observedFrequencies: [frequency, frequency],
        ...jResAnalyzer,
      }),
    );
  }

  return {
    isTwoD: true,
    name: data.source.name,
    expno: data.source.expno,
    zones: newZones,
    experimental: minMax,
  };
}
