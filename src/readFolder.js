'use strinct';
const { writeFileSync } = require('fs');
const os = require('os');

const { gsd, optimizePeaks } = require('ml-gsd');
const { xyExtract } = require('ml-spectra-processing');
const { unparse } = require('papaparse');

const convertSpectra = require('./util/convertSpectra');
const getFolders = require('./util/getFolders');

const sqrtPI = Math.sqrt(Math.PI);

let separator = os.type() === 'Windows_NT' ? '\\' : '/';

let path = '/data2/BIOGUNE';
// let path = 'C:\\users\\alejo\\documents\\BIOGUNE';

let ROI = [
  {
    name: 'tsp',
    range: { from: -0.1, to: 0.1 },
    delta: 0,
    pattern: [1],
    byCandidate: false,
  },
  {
    name: 'eretic',
    range: { from: 14.5, to: 15.5 },
    delta: 15,
    pattern: [1],
    byCandidate: false,
  },
  {
    name: 'glucose',
    delta: 5.23,
    range: { from: 5.15, to: 5.35 },
    pattern: [1, 1.008460237],
    integral: [1],
    jCoupling: [3.74],
    byCandidate: true,
  },
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
    kind: 'pseudoVoigt',
  },
  optimization: {
    kind: 'lm',
    options: {
      maxInterations: 500,
      gradientDifference: 1e-3,
    },
  },
};

let { folders, quantFactorSample } = getFolders(path, { separator });

let result = [];
for (let folder of folders) {
  console.log(`\n\n ${folder}`);
  let spectra = convertSpectra(folder, { separator, xy: true });
  let ereticFactor = quantFactorSample[spectra.filename];
  let {
    data: spectrum,
    deltaX,
    nbPoints,
    observeFrequency: field,
  } = spectra.value.spectra[0];

  if (deltaX < 0) {
    spectrum.x.reverse();
    spectrum.y.reverse();
  }

  //normalize with ereticFactor
  for (let i = 0; i < nbPoints; i++) {
    spectrum.y[i] /= ereticFactor;
  }

  let roiResult = {};
  let first;
  for (let roi of ROI) {
    console.log(`roi: ${roi.name}`);
    let peaks = getOptPeaks(spectrum, {
      gsdOptions,
      optimizationOptions,
      roi: roi.range,
    });

    let bestCandidate, candidates;
    if (roi.byCandidate) {
      // console.log(peaks)
      candidates = lookingForCandidates(peaks, roi, { field });
      candidates.sort((a, b) => b.score - a.score);
      bestCandidate = candidates[0].peaks;
    } else {
      candidates = { peaks: peaks.slice(), score: 1 };
      peaks.sort((a, b) => b.y - a.y);
      bestCandidate = peaks.slice(0, roi.pattern.length);
    }

    let shift =
      bestCandidate.reduce((a, b) => a + b.x, 0) / bestCandidate.length;

    let integration = bestCandidate.reduce((totalIntegration, peak) => {
      let { y, width, mu } = peak;
      totalIntegration += y * width * sqrtPI * (1 - mu + mu * sqrtPI);
      return totalIntegration;
    }, 0);

    roiResult[roi.name] = {
      shift,
      integration,
      candidates,
      peaks: bestCandidate,
    };
  }
  result.push({
    name: spectra.filename,
    rois: roiResult,
  });
}

writeFileSync('peakResult.json', JSON.stringify(result));
writeFileSync('peakResult.csv', unparse(result));

function getOptPeaks(spectrum, options = {}) {
  let { roi, optimizationOptions, gsdOptions } = options;
  let { from, to } = roi;
  let experimental = xyExtract(spectrum, { zones: [{ from, to }] });
  let peaks = gsd(experimental, gsdOptions);
  let optPeaks = optimizePeaks(experimental, peaks, optimizationOptions);
  return optPeaks.map((peak) => {
    let { x, y, width, mu } = peak;
    return { x, y, width, mu };
  });
}
function getCandidates(peaks, jcp, pattern, candidates, options = {}) {
  if (candidates.length === 0) return null;
  if (
    pattern.length === 0 ||
    candidates.some((e) => e.indexs.length === pattern.length)
  ) {
    let { delta, range, nH } = options;
    return candidates.map((cand) => {
      let indexs = cand.indexs;
      let toExport = {
        peaks: indexs.map((index) => peaks[index]),
        score: cand.score,
      };
      return toExport;
    });
  }
  let len = peaks.length;
  let newCandidates = [];
  for (let i = 0; i < candidates.length; i++) {
    let { indexs, score } = candidates[i];
    let index = indexs[indexs.length - 1];
    let iPattern = indexs.length - 1;
    let maxDiff = jcp[iPattern] + 0.001;
    for (let j = index + 1; j < len; j++) {
      let c = Math.abs(peaks[index].x - peaks[j].x);
      let diff = Math.abs(c - jcp[iPattern]) / jcp[iPattern];

      if (c > maxDiff) {
        break;
      }
      if (diff < 0.1) {
        let RIP = pattern[iPattern] / pattern[iPattern + 1];
        let RIC = peaks[index].y / peaks[j].y;

        let diffRI = Math.abs(RIP - RIC) / RIP;

        if (diffRI < 0.08) {
          score += 1 - diffRI;
          newCandidates.push({ indexs: indexs.concat([j]), score });
        }
      } else {
      }
    }
  }
  return getCandidates(peaks, jcp, pattern, newCandidates, options);
}
function lookingForCandidates(peaks, signal, options) {
  let { field } = options;
  let delta = signal.delta;
  let pattern = signal.pattern;
  let jCoupling = signal.jCoupling.map((j) => j / field);
  let candidates = [];
  let candOptions = { delta, range: signal.range, nH: signal.integral };
  peaks.forEach((_peak, pp, arr) => {
    let cand = getCandidates(
      arr,
      jCoupling,
      pattern,
      [{ indexs: [pp], score: 0 }],
      candOptions,
    ); // generate combinations from J
    if (cand !== null) candidates = candidates.concat(cand);
  });
  return candidates;
}
