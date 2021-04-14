const { readFileSync } = require('fs');

const { convertFolder } = require('brukerconverter');
const { IOBuffer } = require('iobuffer');
const walker = require('klaw-sync');

const BINARY = 1;
const TEXT = 2;
const files = {
  ser: BINARY,
  fid: BINARY,
  acqus: TEXT,
  acqu2s: TEXT,
  procs: TEXT,
  proc2s: TEXT,
  '1r': BINARY,
  '1i': BINARY,
  '2rr': BINARY,
};

module.exports = function convertSpectra(folderName, options = {}) {
  let { separator } = options;
  let inFilter = (dir) => {
    let path = dir.path;
    let relativePath = path.substr(
      path.lastIndexOf(separator) + 1,
      path.length,
    );
    return files[relativePath] ? true : false;
  };

  let currFolder = walker(folderName);
  let currFiles = currFolder.filter(inFilter).map((e) => e.path);
  let brukerFiles = {};
  if (folderName.indexOf('pdata') >= 0) {
    brukerFiles.acqus = readFileSync(
      folderName.replace(
        folderName.substr(folderName.lastIndexOf('pdata'), folderName.length),
        'acqus',
      ),
      'utf8',
    );
  }
  for (let j = 0; j < currFiles.length; ++j) {
    let idx = currFiles[j].lastIndexOf(separator);
    let name = currFiles[j].substr(idx + 1);
    if (files[name] === BINARY) {
      brukerFiles[name] = new IOBuffer(readFileSync(currFiles[j]));
    } else {
      brukerFiles[name] = readFileSync(currFiles[j], 'utf8');
    }
  }
  return {
    filename: folderName,
    value: convertFolder(brukerFiles, options),
  };
};
