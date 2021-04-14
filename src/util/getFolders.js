const { join } = require('path');
const walker = require('klaw-sync');
const { readFileSync, existsSync } = require('fs');
const { create: createStream } = require('xml-reader');

module.exports = function getFolders(pathFolder, options = {}) {
  let { separator } = options;

  let folders = [];
  let quantFactorSample = {};
  let directories = walker(pathFolder, { nofile: true });
  for (let i = 0; i < directories.length; i++) {
    let dir = directories[i].path;
console.log(join(dir, '1r'))
    if (!existsSync(join(dir, '1r'))) continue;
//console.log(`pasa first exist`)
    let parts = dir.split(separator);
    if (Number(parts[parts.length - 1]) > 10) continue;
    let quantFactorPath = join(
      parts[0][0] !== dir[0] ? '/':'',
      ...parts.slice(0, parts.length - 2),
      'QuantFactorSample.xml',
    );
//console.log(`path q ${quantFactorPath}`);
    if (!existsSync(quantFactorPath)) continue;
//    console.log('pasa quant exist');
    quantFactorSample[dir] = extractEreticFactor(readFileSync(quantFactorPath, 'utf8'));
    folders.push(dir);
  }
//console.log(folders)
  return { folders, quantFactorSample };
};

function extractEreticFactor(xml) {
  let ereticFactor = 0;
  let reader = createStream({stream: true});
  reader.on('tag:Eretic_Factor', (data) =>  {
    if (data.parent.name === 'Application_Parameter') ereticFactor = data.children[0].value;
  });
  reader.parse(xml);
  return Number(ereticFactor);
}
