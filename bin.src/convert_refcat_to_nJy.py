#!/usr/bin/env python
# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Convert old HTM reference catalogs to use nJy for fluxes, instead of Jy.

This will process all .fits files in the given directory, converting the units
of any `_flux`, `_fluxSigma`, and `_fluxErr` fields from Jy->nJy and making
their units field in the schema `nJy`, then overwriting the files with the
updated values. Also replaces `_fluxSigma` with `_fluxErr` in the schema,
per RFC-333. If no such flux fields are found with units of `''`, `'?'`, or
`'Jy'`, the file is not modified.

Many of our old reference catalogs have no units for their fluxes: we assume
(as the algorithmic code did) that these are all in Jy units.

By default, only print the files and fields that will be modified, but do not
write the output.

If you are processing a large number of files (e.g. ps1_pv3), we recommend
capturing stdout to a log file.
"""
import os.path
import glob

from concurrent.futures import ProcessPoolExecutor
import itertools

import lsst.afw.table


def process_one(filename, write=False):
    """Convert one file in-place from Jy (or no units) to nJy fluxes.
    """
    status = []
    status.append(f"Reading: {filename}")
    catalog = lsst.afw.table.SimpleCatalog.readFits(filename)

    def is_flux_field(name, units):
        unitsCheck = (oldUnits == 'Jy' or oldUnits == '' or oldUnits == '?')
        isFlux = name.endswith('_flux')
        isFluxSigma = name.endswith('_fluxSigma')
        isFluxErr = name.endswith('_fluxErr')
        return (isFlux or isFluxSigma or isFluxErr) and unitsCheck

    # Do not share the AliasMap: for refcats, that gets created when the
    # catalog is read from disk, and should not be propagated in the files.
    mapper = lsst.afw.table.SchemaMapper(catalog.schema, shareAliasMap=False)
    mapper.addMinimalSchema(lsst.afw.table.SimpleTable.makeMinimalSchema())
    mapped_fluxes = []
    for field in catalog.schema:
        oldName = field.field.getName()
        oldUnits = field.field.getUnits()
        if is_flux_field(oldName, oldUnits):
            units = 'nJy'
            # remap Sigma flux fields to Err, since the afw table logic for that mapping
            # seems to rely on a lack of units...
            if oldName.endswith('_fluxSigma'):
                name = oldName.strip('_fluxSigma') + '_fluxErr'
            else:
                name = oldName
            newField = lsst.afw.table.Field[field.dtype](name, field.field.getDoc(), units)
            mapper.addMapping(field.getKey(), newField)
            mapped_fluxes.append((oldName, oldUnits))
        else:
            mapper.addMapping(field.getKey())

    if mapped_fluxes == []:
        status.append("No fluxes to map: not writing.")
        return status

    flux_fields_str = '; '.join("%s, '%s'" % (name, units) for name, units in mapped_fluxes)
    status.append(f"Found flux fields: {flux_fields_str}")

    if write:
        newSchema = mapper.getOutputSchema()
        output = lsst.afw.table.SimpleCatalog(newSchema)
        output.extend(catalog, mapper=mapper)
        for field in output.schema:
            name = field.field.getName()
            if name.endswith('_flux') or name.endswith('_fluxErr'):
                output[field.key] *= 1e9

        output.writeFits(filename)
        status.append(f"Wrote: {filename}")

    return status


def main():
    import argparse

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=CustomFormatter)
    parser.add_argument("path", metavar="path",
                        help="Directory containing the reference catalogs to overwrite."
                        " All files with a `.fits` extension in the directory will be processed.")
    parser.add_argument('-n', '--nprocesses', default=4, type=int,
                        help="Number of processes to use when reading and writing files.")
    parser.add_argument('--write', action="store_true",
                        help="Write the corrected files (default just prints what would have changed).")
    args = parser.parse_args()

    files = glob.glob(os.path.join(args.path, "*.fits"))

    with ProcessPoolExecutor(max_workers=args.nprocesses) as executor:
        mapped = executor.map(process_one, files, itertools.repeat(args.write))
    for result in mapped:
        print('\n'.join(result))


if __name__ == "__main__":
    main()
