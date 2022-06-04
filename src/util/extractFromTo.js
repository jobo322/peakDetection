'use strict';

const { xFindClosestIndex, xSequentialFill } = require('ml-spectra-processing');

module.exports = function extractFromTo(xyData, options = {}) {
  let { from, to, interpolationFct = linearInterpolation, nbPoints } = options;
  const { x, y } = xyData;

  let flip = x[0] > x[1];
  if (flip) {
    x.reverse();
    y.reverse();
  }

  if (from > to) [from, to] = [to, from];

  const fromClosestIndex = xFindClosestIndex(x, from);
  const toClosestIndex = xFindClosestIndex(x, to);

  const newX = xSequentialFill({
    from,
    to,
    size: nbPoints || toClosestIndex - fromClosestIndex,
  });

  const newY = new Float64Array(newX.length);

  let lowerClosestIndex =
    x[fromClosestIndex] > from ? fromClosestIndex - 1 : fromClosestIndex;

  for (let i = 0; i < newX.length; i++) {
    lowerClosestIndex = getLowerClosestIndex(x, newX[i], lowerClosestIndex);
    newY[i] = interpolationFct(x, y, {
      target: newX[i],
      closestLowIndex: lowerClosestIndex,
    });
  }

  if (flip) {
    newX.reverse();
    newY.reverse();
  }

  return { x: newX, y: newY };
};

function linearInterpolation(x, y, options) {
  const { target, closestLowIndex: c } = options;
  return (
    ((y[c + 1] - y[c + 0]) / (x[c + 1] - x[c + 0])) * (target - x[c + 0]) +
    y[c + 0]
  );
}

function getLowerClosestIndex(x, target, lowerClosestIndex) {
  let index =
    xFindClosestIndex(x.slice(lowerClosestIndex), target) + lowerClosestIndex;
  return x[index] > target ? Math.max(index - 1, 0) : index;
}
