'use strict';

const { writeFileSync } = require('fs');
const { join } = require('path');

const { gsd, optimizePeaks } = require('ml-gsd');
const { xyExtract } = require('ml-spectra-processing');
const { SpectrumGenerator } = require('spectrum-generator');
const { solventSuppression } = require('nmr-processing');

const readNMRRSync = require('./util/readNMRRSync');

const path =
  'C:\\Users\\alejo\\Downloads\\Telegram Desktop\\dataset489_469.zip';
const pathToWrite = './';

let ROI = [
  { from: 2.0, to: 2.1 },
  { from: 3.0, to: 3.2 },
];

let gsdOptions = {
  minMaxRatio: 0.01,
  broadRatio: 0.00025,
  smoothY: true,
  realTopDetection: true,
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

const pdata = readNMRRSync(path, { xy: true });

let result = [];
for (let i = 0; i < pdata.length; i++) {
  let filename = pdata[i].filename;
  let spectrum = pdata[i].value.spectra[0].data;
  if (spectrum.x[0] > spectrum.x[1]) {
    spectrum.x = spectrum.x.reverse();
    spectrum.y = spectrum.y.reverse();
  }

  const xyData = align({
    xyData: spectrum,
    referencePeaks: [
      { x: 5.226, y: 1 },
      { x: 5.22, y: 1 },
    ],
    delta: 5.22,
    fromTo: { from: 5.0, to: 5.4 },
  });
  let peakOptimized = { filename, fit: [] };
  for (let roi of ROI) {
    let { from, to } = roi;
    let experimental = xyExtract(xyData, { zones: [{ from, to }] });
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
  result.push(peakOptimized);
}

writeFileSync(join(pathToWrite, 'fittingResult.json'), JSON.stringify(result));

function align(input) {
  const { xyData, referencePeaks, delta, fromTo } = input;

  const ROI = xyExtract(xyData, { zones: [fromTo] });
  const peaks = gsd(ROI, gsdOptions);

  const marketPeaks = solventSuppression(
    peaks,
    {
      delta,
      peaks: referencePeaks,
    },
    { marketPeaks: true },
  );

  const glucosePeaks = marketPeaks.filter((peak) => peak.kind === 'solvent');
  const shift =
    delta - glucosePeaks.reduce((a, b) => a + b.x, 0) / glucosePeaks.length;

  xyData.x.forEach((e, i, arr) => (arr[i] += shift));

  return xyData;
}
