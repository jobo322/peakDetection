'use strict';

const { resolve } = require('path');

const { convertFolder } = require('brukerconverter');
const { IOBuffer } = require('iobuffer');
const { convert } = require('jcampconverter');
const zipper = require('zip-local');

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

module.exports = function readNMRRSync(path, options = {}) {
  let zipData = zipper.sync.unzip(resolve(path)).memory();
  let zipFiles = zipData.unzipped_file;
  let folders = getSpectraFolders(zipFiles);
  let spectra = convertSpectraSync(folders.brukerFolders, zipFiles, options);
  let jcamps = processJcamp(folders.jcampFolders, zipFiles);
  return spectra.concat(jcamps);
};

function convertSpectraSync(folders, zip, options = {}) {
  let spectra = new Array(folders.length);
  let inFilter = (relativePath, file) => {
    return files[relativePath] ? true : false;
  };
  for (let i = 0; i < folders.length; ++i) {
    let len = folders[i].name.length;
    let folderName = folders[i].name;
    folderName = folderName.substr(0, folderName.lastIndexOf('/') + 1);
    let currFolder = zip.folder(folderName);
    let currFiles = currFolder.filter(inFilter);
    let brukerFiles = {};
    if (folderName.indexOf('pdata') >= 0) {
      brukerFiles.acqus = zip
        .file(folderName.replace(/pdata\/[0-9]\//, 'acqus'))
        .asText();
    }
    for (let j = 0; j < currFiles.length; ++j) {
      let idx = currFiles[j].name.lastIndexOf('/');
      let name = currFiles[j].name.substr(idx + 1);
      if (files[name] === BINARY) {
        brukerFiles[name] = new IOBuffer(currFiles[j].asArrayBuffer());
      } else {
        brukerFiles[name] = currFiles[j].asText();
      }
    }
    spectra[i] = {
      filename: folderName,
      value: convertFolder(brukerFiles, options),
    };
  }
  return spectra;
}

function getSpectraFolders(zipFiles) {
  // Folders should contain jcamp too
  let brukerFolders = zipFiles.filter((relativePath) => {
    if (relativePath.match('__MACOSX')) return false;
    if (relativePath.endsWith('1r')) {
      return true;
    }
    return false;
  });
  let jcampFolders = zipFiles.filter((relativePath) => {
    if (relativePath.match('__MACOSX')) return false;
    if (relativePath.endsWith('dx') || relativePath.endsWith('jcamp')) {
      return true;
    }
    return false;
  });
  return { jcampFolders, brukerFolders };
}

function processJcamp(folders, zipFiles, options) {
  let spectra = new Array(folders.length);
  for (let i = 0; i < folders.length; ++i) {
    let name = folders[i].name;
    let jcamp = zipFiles.file(name).asText();
    let value = convert(jcamp, {
      keepSpectra: true,
      keepRecordsRegExp: /^.+$/,
      xy: true,
    });
    spectra[i] = { filename: name, value };
  }
  return spectra;
}
