'use strict';

const { readFileSync, existsSync } = require('fs');
const { join } = require('path');

const walker = require('klaw-sync');
const { create: createStream } = require('xml-reader');

module.exports = function getFolders(pathFolder, options = {}) {
  let { separator } = options;
  let folders = [];
  let quantFactorSample = {};
  let first = false;
  let directories = walker(pathFolder, { nofile: true });
  for (let i = 0; i < directories.length; i++) {
    if (first) continue;
    let dir = directories[i].path;
    if (!existsSync(join(dir, '1r'))) continue;
    let parts = dir.split(separator);
    if (Number(parts[parts.length - 1]) > 10) continue;
    let quantFactorPath = join(
      parts[0][0] !== dir[0] ? '/' : '',
      ...parts.slice(0, parts.length - 2),
      'QuantFactorSample.xml',
    );
    if (!existsSync(quantFactorPath)) continue;
    quantFactorSample[dir] = extractEreticFactor(
      readFileSync(quantFactorPath, 'utf8'),
    );
    folders.push(dir);
    // first = true;
  }
  return { folders, quantFactorSample };
};

function extractEreticFactor(xml) {
  let result = {};
  let reader = createStream({ stream: true });
  reader.on('tag:Eretic_Factor', (data) => {
    if (data.parent.name === 'Application_Parameter') {
      result.ereticFactor = Number(data.children[0].value);
    }
  });
  reader.on('tag:Experiment_Description', (data) => {
    if (data.parent.name === 'Application_Parameter') {
      result.experimentDescription = toJSON(data, {});
    }
  });
  reader.parse(xml);

  let { experimentDescription, ereticFactor } = result;
  return { ereticFactor, ...experimentDescription };
}

function toJSON(tag, result) {
  if (!tag.name && tag.children.length === 0) return tag.value;
  for (let child of tag.children) {
    let { children, value, name } = child;
    if (name.length > 0) {
      result[name] = toJSON(child, {});
      continue;
    }
    if (name.length === 0 && children.length === 0) return value;
  }
  return result;
}
