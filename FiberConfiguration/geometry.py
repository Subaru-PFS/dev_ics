import argparse
import logging
import os
import pathlib
import sys
import warnings
from datetime import datetime

import h5py
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

# Load site-specific modules only after standard environment packages.
sys.path.append(f'{os.environ["HOME"]}/pythonmodules')
sys.path.append(f'{os.environ["HOME"]}/pythonmodules/selenium')

from ics.cobraCharmer.cobraCoach import visDianosticPlot_TK as vdp


DEFAULT_XML = (
    "/software/mhs/products/Linux64/"
    "pfs_instdata/1.7.36/data/pfi/modules/ALL/ALL.xml"
)

COMMON_GEOMETRY_DATASETS = (
    "newCnts_circ",
    "newRs_circ",
    "CircRes",
)

THETA_GEOMETRY_DATASETS = (
    *COMMON_GEOMETRY_DATASETS,
    "newCnts_ellip",
    "newRs_ellip",
    "PA_ellip",
    "EllipRes",
)


def runDir2Date(runDir):
    """Return the date encoded by a yyyymmdd_numeric run directory."""

    try:
        dateText, runNumber = str(runDir).split("_", maxsplit=1)

        if not runNumber or not runNumber.isdigit():
            raise ValueError

        return datetime.strptime(dateText, "%Y%m%d")
    except ValueError as exc:
        raise ValueError(
            f"Invalid runDir {runDir!r}; expected yyyymmdd followed "
            "by an underscore and a nonempty numeric run number."
        ) from exc


def normalizeEL(value):
    """Return a finite numeric elevation, preserving integral values."""

    if isinstance(value, bytes):
        value = value.decode("utf-8")

    try:
        elevation = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid EL value: {value!r}") from exc

    if not np.isfinite(elevation):
        raise ValueError(f"Invalid EL value: {value!r}")

    return int(elevation) if elevation.is_integer() else elevation


def askYesNo(prompt):
    while True:
        answer = input(prompt).strip().lower()

        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False

        print("Please type yes or no.")


def visit2runDir(
    visit,
    arm,
    EL,
    *,
    visitRunDirList="visit_runDir.csv",
    update=True,
):
    """Resolve a visit to a runDir and maintain its CSV metadata."""

    arm = str(arm).lower()
    if arm not in {"theta", "phi"}:
        raise ValueError("arm must be either 'theta' or 'phi'.")

    EL = normalizeEL(EL)
    csvPath = pathlib.Path(visitRunDirList)

    try:
        table = pd.read_csv(csvPath)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Could not find {csvPath}.") from exc

    required = {"visit", "runDir", "arm", "EL"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            f"{csvPath} is missing columns: {sorted(missing)}"
        )

    visitMask = table["visit"].eq(visit)
    matches = table.loc[visitMask]

    if len(matches) == 1:
        runDir = str(matches.iloc[0]["runDir"])
    else:
        if len(matches) > 1:
            warnings.warn(
                f"Multiple runDir entries exist for visit {visit}; "
                "resolving it again with vdp.findRunDir().",
                RuntimeWarning,
            )

        runDir = vdp.findRunDir(visit)
        if runDir is None:
            raise FileNotFoundError(
                f"No runDir could be found for visit {visit}."
            )
        runDir = str(runDir)

    runDir2Date(runDir)
    selectedMask = visitMask & table["runDir"].astype(str).eq(runDir)
    changed = False

    if not selectedMask.any():
        runDate = runDir2Date(runDir)
        newRow = {column: np.nan for column in table.columns}
        newRow.update(
            visit=visit,
            runDir=runDir,
            arm=arm,
            EL=EL,
        )
        if "year" in table.columns:
            newRow["year"] = runDate.year
        if "month" in table.columns:
            newRow["month"] = runDate.month

        table = pd.concat(
            [table, pd.DataFrame([newRow])],
            ignore_index=True,
        )
        selectedMask = (
            table["visit"].eq(visit)
            & table["runDir"].astype(str).eq(runDir)
        )
        changed = True

    selected = table.loc[selectedMask]
    oldArm = str(selected.iloc[0]["arm"]).lower()
    oldEL = normalizeEL(selected.iloc[0]["EL"])

    shouldUpdateArm = oldArm != arm and (
        update
        or askYesNo(
            f"CSV arm is {oldArm!r}, requested arm is {arm!r}. "
            "Update it? [y/n] "
        )
    )
    shouldUpdateEL = oldEL != EL and (
        update
        or askYesNo(
            f"CSV EL is {oldEL!r}, requested EL is {EL!r}. "
            "Update it? [y/n] "
        )
    )

    if shouldUpdateArm:
        table.loc[selectedMask, "arm"] = arm
        changed = True
    if shouldUpdateEL:
        table.loc[selectedMask, "EL"] = EL
        changed = True

    duplicateMask = visitMask & ~table["runDir"].astype(str).eq(runDir)
    if duplicateMask.any():
        removeDuplicates = update or askYesNo(
            f"Remove other runDir rows for visit {visit}? [y/n] "
        )
        if removeDuplicates:
            table = table.loc[~duplicateMask].copy()
            changed = True

    if changed:
        table.to_csv(csvPath, index=False)

    return runDir


def runDir2Visit(runDir, *, visitRunDirList="visit_runDir.csv"):
    """Resolve exactly one visit associated with a runDir in the CSV."""

    csvPath = pathlib.Path(visitRunDirList)
    try:
        table = pd.read_csv(csvPath)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Could not find {csvPath}.") from exc

    required = {"visit", "runDir"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            f"{csvPath} is missing columns: {sorted(missing)}"
        )

    matches = table.loc[
        table["runDir"].astype(str).eq(str(runDir)),
        "visit",
    ].dropna()
    visits = pd.unique(matches)

    if len(visits) == 0:
        raise LookupError(
            f"No visit is recorded for runDir {runDir!r} in {csvPath}."
        )
    if len(visits) > 1:
        raise ValueError(
            f"Multiple visits are recorded for runDir {runDir!r} "
            f"in {csvPath}: {visits.tolist()}"
        )

    try:
        return int(visits[0])
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid visit {visits[0]!r} for runDir {runDir!r}."
        ) from exc


def visitList2runDirList(
    visitsFile,
    *,
    visitRunDirList="visit_runDir.csv",
    update=True,
    errorHandler=None,
):
    """Read `runNumber visit arm EL rot` rows and resolve runDirs."""

    visitTable = np.loadtxt(visitsFile, dtype=str, ndmin=2)
    if visitTable.shape[1] != 5:
        raise ValueError(
            f"{visitsFile} must contain: "
            "runNumber, visit, arm, EL, and rot."
        )

    rowCount = len(visitTable)
    runNumbers = np.full(rowCount, None, dtype=object)
    visits = np.full(rowCount, None, dtype=object)
    arms = np.full(rowCount, None, dtype=object)
    elevations = np.full(rowCount, np.nan, dtype=float)
    rotations = np.full(rowCount, np.nan, dtype=float)
    runDirs = np.full(rowCount, None, dtype=object)

    for index, row in enumerate(visitTable):
        try:
            runNumber = int(row[0])
            visit = int(row[1])
            arm = str(row[2]).lower()
            if arm not in {"theta", "phi"}:
                raise ValueError(
                    f"arm must be theta or phi, not {row[2]!r}"
                )
            elevation = normalizeEL(row[3])
            rotation = normalizeEL(row[4])

            runNumbers[index] = runNumber
            visits[index] = visit
            arms[index] = arm
            elevations[index] = elevation
            rotations[index] = rotation
            runDirs[index] = visit2runDir(
                visit,
                arm,
                elevation,
                visitRunDirList=visitRunDirList,
                update=update,
            )
        except Exception as exc:
            item = f"row {index + 1} ({' '.join(row)})"
            if errorHandler is None:
                warnings.warn(
                    f"Could not process {item}: {exc}",
                    RuntimeWarning,
                )
            else:
                errorHandler(item, exc)

    return runNumbers, runDirs, visits, arms, elevations, rotations


def circleFit(
    runDir,
    arm,
    EL,
    *,
    runNumber=None,
    rot=None,
    xml=DEFAULT_XML,
    thetaHdf5="geometry/thetaGeometry.hdf5",
    phiHdf5="geometry/phiGeometry.hdf5",
    update=True,
):
    """Fit cobra circle geometry and theta ellipse geometry."""

    arm = str(arm).lower()
    if arm == "theta":
        hdf5File = thetaHdf5
    elif arm == "phi":
        hdf5File = phiHdf5
    else:
        raise ValueError("arm must be either 'theta' or 'phi'.")

    runDir2Date(runDir)
    EL = normalizeEL(EL)
    if (runNumber is None) != (rot is None):
        raise ValueError("runNumber and rot must be supplied together.")
    if runNumber is not None:
        runNumber = int(runNumber)
        rot = normalizeEL(rot)

    try:
        vis = vdp.VisDianosticPlot(
            runDir=runDir,
            xml=pathlib.Path(xml),
            arm=arm,
            datatype="MM",
        )
        vis.fw
    except Exception as exc:
        raise AttributeError(
            f"Cannot import {runDir} Cobra data."
        ) from exc

    cobraCount = vis.fw.shape[0]
    nanComplex = np.nan + 1j * np.nan
    newCenters = np.full((cobraCount, 2), nanComplex, dtype=complex)
    newRadii = np.full((cobraCount, 2), nanComplex, dtype=complex)
    positionAngles = np.full(cobraCount, np.nan, dtype=float)
    residuals = np.full((cobraCount, 2), np.nan, dtype=float)

    previousLoggingLevel = logging.root.manager.disable
    logging.disable(logging.CRITICAL)

    try:
        for index in tqdm(
            range(cobraCount),
            desc=f"Performing fits for {runDir}",
        ):
            fw = vis.fw[index, 0, :]
            rv = vis.rv[index, 0, :]
            data = np.concatenate([fw, rv])

            if np.std(fw) < 1 or np.std(rv) < 1:
                continue

            if arm == "theta":
                estimatedCenter = vis.calibModel.centers[index]
                if not (
                    np.isfinite(estimatedCenter.real)
                    and np.isfinite(estimatedCenter.imag)
                ):
                    estimatedCenter = np.median(data)
            else:
                estimatedCenter = np.median(data)

            circle = vdp.ransac_circle(data, estimatedCenter)
            newCenters[index, 0], newRadii[index, 0] = circle.fit(
                ransacTolerance=0.1
            )
            residuals[index, 0] = circle.best_res

            if arm == "theta":
                ellipse = vdp.ransac_ellipse(data, estimatedCenter)
                (
                    newCenters[index, 1],
                    newRadii[index, 1],
                    positionAngles[index],
                ) = ellipse.fit(ransacTolerance=0.1)
                residuals[index, 1] = ellipse.best_res
    finally:
        logging.disable(previousLoggingLevel)

    outputPath = pathlib.Path(hdf5File)
    outputPath.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(outputPath, "a") as hdf5:
        groupName = str(runDir)
        if groupName in hdf5:
            replace = update or askYesNo(
                f"Update {runDir} in {outputPath}? [y/n] "
            )
            if not replace:
                return runDir
            del hdf5[groupName]

        group = hdf5.create_group(groupName)
        group.attrs["runDir"] = str(runDir)
        group.attrs["EL"] = EL
        if runNumber is not None:
            group.attrs["runNumber"] = runNumber
            group.attrs["rot"] = rot

        cobraIdx = np.sort(
            np.concatenate([vis.goodIdx, vis.badIdx])
        )
        if len(cobraIdx) != cobraCount:
            cobraIdx = np.arange(cobraCount)

        group.create_dataset("cobraIdx", data=cobraIdx)
        group.create_dataset(
            "newCnts_circ",
            data=np.column_stack(
                [newCenters[:, 0].real, newCenters[:, 0].imag]
            ),
        )
        group.create_dataset("newRs_circ", data=newRadii[:, 0].real)
        group.create_dataset("CircRes", data=residuals[:, 0])

        if arm == "theta":
            group.create_dataset(
                "newCnts_ellip",
                data=np.column_stack(
                    [newCenters[:, 1].real, newCenters[:, 1].imag]
                ),
            )
            group.create_dataset(
                "newRs_ellip",
                data=np.column_stack(
                    [newRadii[:, 1].real, newRadii[:, 1].imag]
                ),
            )
            group.create_dataset("PA_ellip", data=positionAngles)
            group.create_dataset("EllipRes", data=residuals[:, 1])

    return runDir


def ensureHdf5File(hdf5File):
    path = pathlib.Path(hdf5File)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "a"):
        pass


def circularNanMedian(values, axis=0, period=np.pi):
    """NaN-aware sample circular median with the requested period."""

    values = np.asarray(values, dtype=float)
    moved = np.moveaxis(values, axis, 0)
    flattened = moved.reshape(moved.shape[0], -1)
    result = np.full(flattened.shape[1], np.nan)

    for columnIndex in range(flattened.shape[1]):
        valid = flattened[:, columnIndex]
        valid = valid[np.isfinite(valid)]
        if valid.size == 0:
            continue

        wrapped = np.mod(valid, period)
        candidates = np.unique(wrapped)
        differences = (
            candidates[:, None] - wrapped[None, :] + period / 2
        ) % period - period / 2
        result[columnIndex] = candidates[
            np.argmin(np.sum(np.abs(differences), axis=1))
        ]

    return result.reshape(moved.shape[1:])


def collectRunsByRunELRot(hdf5File, datasetNames, arm):
    groupedRuns = {}

    with h5py.File(hdf5File, "r") as hdf5:
        for runDir, group in hdf5.items():
            if not isinstance(group, h5py.Group):
                continue
            missingAttrs = {
                "runNumber", "EL", "rot"
            } - set(group.attrs.keys())
            if missingAttrs:
                warnings.warn(
                    f"Skipping {hdf5File}/{runDir}: missing attributes "
                    f"{sorted(missingAttrs)}.",
                    RuntimeWarning,
                )
                continue

            runNumber = int(group.attrs["runNumber"])
            EL = normalizeEL(group.attrs["EL"])
            rot = normalizeEL(group.attrs["rot"])
            key = (runNumber, EL, rot, arm)
            required = {"cobraIdx", *datasetNames}
            missing = required - set(group.keys())
            if missing:
                raise KeyError(
                    f"{hdf5File}/{runDir} is missing: {sorted(missing)}"
                )

            cobraIdx = np.asarray(group["cobraIdx"][:])
            if cobraIdx.ndim != 1:
                raise ValueError(
                    f"{hdf5File}/{runDir}/cobraIdx must be 1-D."
                )
            if len(np.unique(cobraIdx)) != len(cobraIdx):
                raise ValueError(
                    f"{hdf5File}/{runDir}/cobraIdx has duplicates."
                )

            data = {}
            for datasetName in datasetNames:
                values = np.asarray(group[datasetName][:], dtype=float)
                if values.ndim == 0 or values.shape[0] != len(cobraIdx):
                    raise ValueError(
                        f"{hdf5File}/{runDir}/{datasetName} does not "
                        "match cobraIdx."
                    )
                data[datasetName] = values

            groupedRuns.setdefault(key, []).append(
                {
                    "runDir": runDir,
                    "date": runDir2Date(runDir),
                    "cobraIdx": cobraIdx,
                    "data": data,
                }
            )

    return groupedRuns


def calculateGeometryMedian(runs, datasetNames):
    allCobraIdx = np.unique(
        np.concatenate([run["cobraIdx"] for run in runs])
    )
    medianData = {}

    for datasetName in datasetNames:
        trailingShapes = {
            run["data"][datasetName].shape[1:] for run in runs
        }
        if len(trailingShapes) != 1:
            raise ValueError(
                f"Inconsistent shapes for {datasetName}: {trailingShapes}"
            )

        trailingShape = trailingShapes.pop()
        stack = np.full(
            (len(runs), len(allCobraIdx), *trailingShape),
            np.nan,
        )

        for runIndex, run in enumerate(runs):
            positions = np.searchsorted(allCobraIdx, run["cobraIdx"])
            stack[runIndex, positions, ...] = run["data"][datasetName]

        if datasetName == "PA_ellip":
            medianData[datasetName] = circularNanMedian(
                stack, axis=0, period=np.pi
            )
        else:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="All-NaN slice encountered",
                    category=RuntimeWarning,
                )
                medianData[datasetName] = np.nanmedian(stack, axis=0)

    return allCobraIdx, medianData


def mergeGeometryByRunELRotArm(
    thetaFile,
    phiFile,
    outputFile,
    *,
    overwrite=True,
):
    ensureHdf5File(thetaFile)
    ensureHdf5File(phiFile)

    thetaRuns = collectRunsByRunELRot(
        thetaFile, THETA_GEOMETRY_DATASETS, "theta"
    )
    phiRuns = collectRunsByRunELRot(
        phiFile, COMMON_GEOMETRY_DATASETS, "phi"
    )

    outputPath = pathlib.Path(outputFile)
    outputPath.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if overwrite else "x"

    with h5py.File(outputPath, mode) as output:
        allRuns = thetaRuns | phiRuns
        keys = sorted(
            allRuns,
            key=lambda key: (
                int(key[0]), float(key[1]), float(key[2]), key[3]
            ),
        )

        runDates = {}
        for (runNumber, _EL, _rot, _arm), runs in allRuns.items():
            runDates.setdefault(runNumber, []).extend(
                run["date"] for run in runs
            )

        for runNumber, dates in runDates.items():
            startDate = min(dates)
            endDate = max(dates)
            middleDate = startDate + (endDate - startDate) / 2
            runGroup = output.require_group(f"run{runNumber}")
            runGroup.attrs["runNumber"] = runNumber
            runGroup.attrs["date"] = middleDate.strftime("%Y%m")
            runGroup.attrs["startDate"] = startDate.strftime("%Y%m%d")
            runGroup.attrs["endDate"] = endDate.strftime("%Y%m%d")

        for runNumber, EL, rot, arm in keys:
            runGroup = output[f"run{runNumber}"]
            elGroup = runGroup.require_group(f"EL_{EL}")
            elGroup.attrs["EL"] = EL
            rotGroup = elGroup.require_group(f"rot_{rot}")
            rotGroup.attrs["rot"] = rot

            groupedRuns = thetaRuns if arm == "theta" else phiRuns
            datasetNames = (
                THETA_GEOMETRY_DATASETS
                if arm == "theta"
                else COMMON_GEOMETRY_DATASETS
            )
            runs = groupedRuns[(runNumber, EL, rot, arm)]
            cobraIdx, medians = calculateGeometryMedian(runs, datasetNames)
            armGroup = rotGroup.create_group(arm)
            armGroup.attrs["arm"] = arm
            armGroup.attrs["numberOfRuns"] = len(runs)
            armGroup.attrs["runDirs"] = np.asarray(
                [run["runDir"] for run in runs],
                dtype=h5py.string_dtype("utf-8"),
            )
            armGroup.create_dataset("cobraIdx", data=cobraIdx)
            for datasetName, values in medians.items():
                armGroup.create_dataset(datasetName, data=values)

    return outputPath


def calCombinedGeometry(geometryFile):
    with h5py.File(geometryFile, "a") as geometry:
        for runName in geometry.keys():
            for elName in geometry[runName].keys():
                elGroup = geometry[runName][elName]
                if not isinstance(elGroup, h5py.Group):
                    continue
                for rotName in elGroup.keys():
                    rotGroup = elGroup[rotName]
                    if "theta" not in rotGroup or "phi" not in rotGroup:
                        warnings.warn(
                            f"Skipping {runName}/{elName}/{rotName}: "
                            "both arms are required.",
                            RuntimeWarning,
                        )
                        continue

                    theta = rotGroup["theta"]
                    phi = rotGroup["phi"]
                    cobraIdx, thetaPositions, phiPositions = np.intersect1d(
                        theta["cobraIdx"][:],
                        phi["cobraIdx"][:],
                        return_indices=True,
                    )
                    if cobraIdx.size == 0:
                        continue

                    thetaCenters = theta["newCnts_circ"][:][thetaPositions]
                    phiCenters = phi["newCnts_circ"][:][phiPositions]
                    thetaRadii = theta["newRs_circ"][:][thetaPositions]
                    phiRadii = phi["newRs_circ"][:][phiPositions]

                    l1 = np.linalg.norm(thetaCenters - phiCenters, axis=1)
                    l2 = phiRadii
                    denominator = 2 * thetaRadii * l1
                    cosine = np.full(l1.shape, np.nan)
                    valid = (
                        np.isfinite(thetaRadii)
                        & np.isfinite(l1)
                        & np.isfinite(l2)
                        & (denominator != 0)
                        & (l2 <= thetaRadii + l1)
                        & (l2 >= np.abs(thetaRadii - l1))
                    )
                    cosine[valid] = (
                        thetaRadii[valid] ** 2
                        + l1[valid] ** 2
                        - l2[valid] ** 2
                    ) / denominator[valid]
                    cosine[valid] = np.clip(cosine[valid], -1, 1)

                    datasets = {
                        "cobraIdx": cobraIdx,
                        "l1": l1,
                        "l2": l2,
                        "opening_angle": np.arccos(cosine),
                    }
                    for name, values in datasets.items():
                        if name in rotGroup:
                            del rotGroup[name]
                        rotGroup.create_dataset(name, data=values)


def pos2angle(x, y):
    """Return a polar angle in the interval [0, 2 pi)."""

    return np.mod(np.arctan2(y, x), 2 * np.pi)


def pos2phi(thetaCntx, thetaCnty, phiCntx, phiCnty, x, y):
    """Return the signed phi angle in the interval [-pi/2, 3pi/2)."""

    reference = np.asarray(
        [thetaCntx - phiCntx, thetaCnty - phiCnty],
        dtype=float,
    )
    position = np.asarray([x - phiCntx, y - phiCnty], dtype=float)

    if not np.all(np.isfinite(reference)) or not np.all(
        np.isfinite(position)
    ):
        return np.nan
    if np.linalg.norm(reference) == 0 or np.linalg.norm(position) == 0:
        return np.nan

    cross = reference[0] * position[1] - reference[1] * position[0]
    angle = np.arctan2(cross, np.dot(reference, position))
    return (angle + 0.5 * np.pi) % (2 * np.pi) - 0.5 * np.pi


def _mergedGeometryGroup(geometry, runNumber, EL, rot):
    """Find the merged run/EL/rotation geometry group."""

    runName = f"run{int(runNumber)}"
    if runName not in geometry:
        raise KeyError(f"No merged geometry is available for {runName}.")
    requestedEL = normalizeEL(EL)
    requestedRot = normalizeEL(rot)
    for elName, group in geometry[runName].items():
        if not isinstance(group, h5py.Group) or "EL" not in group.attrs:
            continue
        if normalizeEL(group.attrs["EL"]) == requestedEL:
            for rotName, rotGroup in group.items():
                if (
                    isinstance(rotGroup, h5py.Group)
                    and "rot" in rotGroup.attrs
                    and normalizeEL(rotGroup.attrs["rot"]) == requestedRot
                ):
                    return rotGroup

    raise KeyError(
        f"No merged geometry is available for {runName}, "
        f"EL={requestedEL}, rot={requestedRot}."
    )


def _cobraPositions(cobraIdx):
    """Map each cobra identifier to its row in an HDF5 dataset."""

    return {int(cobra): index for index, cobra in enumerate(cobraIdx)}


def _onCircle(point, center, radius, tolerance):
    """Test a measured point against a fitted circle."""

    values = np.asarray(
        [point.real, point.imag, center[0], center[1], radius],
        dtype=float,
    )
    if not np.all(np.isfinite(values)):
        return False

    radiusResidual = (
        (point.real - center[0]) ** 2
        + (point.imag - center[1]) ** 2
        - radius**2
    )
    return abs(radiusResidual) <= tolerance


def _replaceDataset(group, name, values, update):
    if name in group:
        replace = update or askYesNo(
            f"Update {group.name}/{name}? [y/n] "
        )
        if not replace:
            return
        del group[name]
    group.create_dataset(name, data=values)


def calculateMovementLimits(
    thetaFile,
    phiFile,
    geometryFile,
    *,
    xml=DEFAULT_XML,
    outlierTolerance=1.0,
    update=True,
    errorHandler=None,
):
    """Calculate run-level CCW, CW, Phi_in, and Phi_out in degrees.

    Circle centers, radii, and opening angles come from geometry merged
    by run number, elevation, rotation, and arm. Results are written into
    corresponding runDir groups in thetaFile and phiFile.
    """

    if outlierTolerance < 0:
        raise ValueError("outlierTolerance must be nonnegative.")

    ensureHdf5File(thetaFile)
    ensureHdf5File(phiFile)
    geometryPath = pathlib.Path(geometryFile)
    if not geometryPath.exists():
        raise FileNotFoundError(
            f"Merged geometry file {geometryPath} does not exist."
        )

    completed = {"theta": [], "phi": []}
    with h5py.File(geometryPath, "r") as geometry:
        for arm, sourceFile, outputNames in (
            ("theta", thetaFile, ("CCW", "CW")),
            ("phi", phiFile, ("Phi_in", "Phi_out")),
        ):
            with h5py.File(sourceFile, "a") as source:
                for runDir, runGroup in source.items():
                    if not isinstance(runGroup, h5py.Group):
                        continue

                    try:
                        requiredAttrs = {"runNumber", "EL", "rot"}
                        missingAttrs = requiredAttrs - set(runGroup.attrs.keys())
                        if missingAttrs:
                            warnings.warn(
                                f"Skipping {sourceFile}/{runDir}: missing "
                                f"attributes {sorted(missingAttrs)}.",
                                RuntimeWarning,
                            )
                            continue
                        if "cobraIdx" not in runGroup:
                            raise KeyError("the run group has no cobraIdx dataset")

                        merged = _mergedGeometryGroup(
                            geometry,
                            runGroup.attrs["runNumber"],
                            runGroup.attrs["EL"],
                            runGroup.attrs["rot"],
                        )
                        required = {"theta", "phi", "cobraIdx", "opening_angle"}
                        missing = required - set(merged.keys())
                        if missing:
                            raise KeyError(
                                "merged geometry is missing "
                                f"{sorted(missing)}"
                            )

                        theta = merged["theta"]
                        phi = merged["phi"]
                        thetaRows = _cobraPositions(theta["cobraIdx"][:])
                        phiRows = _cobraPositions(phi["cobraIdx"][:])
                        combinedRows = _cobraPositions(merged["cobraIdx"][:])

                        thetaCenters = theta["newCnts_circ"][:]
                        thetaRadii = theta["newRs_circ"][:]
                        phiCenters = phi["newCnts_circ"][:]
                        phiRadii = phi["newRs_circ"][:]
                        openingAngles = merged["opening_angle"][:]

                        vis = vdp.VisDianosticPlot(
                            runDir=runDir,
                            xml=pathlib.Path(xml),
                            arm=arm,
                            datatype="MM",
                        )
                        cobraIdx = np.asarray(runGroup["cobraIdx"][:])
                        first = np.full(len(cobraIdx), np.nan)
                        second = np.full(len(cobraIdx), np.nan)

                        for row, cobra in enumerate(cobraIdx):
                            cobra = int(cobra)
                            if row >= vis.fw.shape[0] or row >= vis.rv.shape[0]:
                                continue

                            fw = vis.fw[row, 0, 0]
                            rv = vis.rv[row, 0, 0]

                            if arm == "theta":
                                if cobra not in thetaRows or cobra not in combinedRows:
                                    continue
                                thetaRow = thetaRows[cobra]
                                combinedRow = combinedRows[cobra]
                                center = thetaCenters[thetaRow]
                                radius = thetaRadii[thetaRow]
                                opening = openingAngles[combinedRow]
                                if not np.isfinite(opening):
                                    continue

                                if _onCircle(fw, center, radius, outlierTolerance):
                                    first[row] = np.degrees(
                                        (pos2angle(
                                            fw.real - center[0],
                                            fw.imag - center[1],
                                        ) + opening) % (2 * np.pi)
                                    )
                                if _onCircle(rv, center, radius, outlierTolerance):
                                    second[row] = np.degrees(
                                        (pos2angle(
                                            rv.real - center[0],
                                            rv.imag - center[1],
                                        ) + opening) % (2 * np.pi)
                                    )
                            else:
                                if cobra not in thetaRows or cobra not in phiRows:
                                    continue
                                thetaCenter = thetaCenters[thetaRows[cobra]]
                                phiRow = phiRows[cobra]
                                phiCenter = phiCenters[phiRow]
                                phiRadius = phiRadii[phiRow]

                                if _onCircle(
                                    fw, phiCenter, phiRadius, outlierTolerance
                                ):
                                    first[row] = np.degrees(
                                        pos2phi(
                                            thetaCenter[0], thetaCenter[1],
                                            phiCenter[0], phiCenter[1],
                                            fw.real, fw.imag,
                                        )
                                    )
                                if _onCircle(
                                    rv, phiCenter, phiRadius, outlierTolerance
                                ):
                                    second[row] = np.degrees(
                                        pos2phi(
                                            thetaCenter[0], thetaCenter[1],
                                            phiCenter[0], phiCenter[1],
                                            rv.real, rv.imag,
                                        )
                                    )

                        _replaceDataset(
                            runGroup, outputNames[0], first, update
                        )
                        _replaceDataset(
                            runGroup, outputNames[1], second, update
                        )
                        completed[arm].append(str(runDir))
                    except Exception as exc:
                        item = (
                            f"{sourceFile}/{runDir} while calculating "
                            "movement limits"
                        )
                        if errorHandler is None:
                            warnings.warn(
                                f"Skipping {item}: {exc}",
                                RuntimeWarning,
                            )
                        else:
                            errorHandler(item, exc)

    return completed


class GeometryMeasurement:
    def __init__(
        self,
        visit=None,
        arm=None,
        EL=None,
        runNumber=None,
        rot=None,
        runDir=None,
        *,
        xml=DEFAULT_XML,
        visitRunDirList="visit_runDir.csv",
        thetaHdf5="geometry/thetaGeometry.hdf5",
        phiHdf5="geometry/phiGeometry.hdf5",
        geometryHdf5="geometry/Geometry.hdf5",
        update=True,
        overwriteCombined=True,
        batchErrorPolicy="warn",
    ):
        if visit is not None and runDir is not None:
            raise ValueError("Specify either visit or runDir, not both.")

        self.visit = visit
        self.arm = str(arm).lower() if arm is not None else None
        self.EL = EL
        self.runNumber = runNumber
        self.rot = rot
        self.runDir = runDir
        self.xml = pathlib.Path(xml)
        self.visitRunDirList = pathlib.Path(visitRunDirList)
        self.thetaHdf5 = pathlib.Path(thetaHdf5)
        self.phiHdf5 = pathlib.Path(phiHdf5)
        self.geometryHdf5 = pathlib.Path(geometryHdf5)
        self.update = update
        self.overwriteCombined = overwriteCombined
        if batchErrorPolicy not in {"warn", "prompt"}:
            raise ValueError(
                "batchErrorPolicy must be 'warn' or 'prompt'."
            )
        self.batchErrorPolicy = batchErrorPolicy
        self.errors = []
        self.runDirs = None
        self.runNumbers = None
        self.visits = None
        self.arms = None
        self.ELs = None
        self.rots = None

    def _validateSingleInput(self):
        if self.arm not in {"theta", "phi"}:
            raise ValueError(
                "arm must be theta or phi for a single visit/runDir."
            )
        if self.EL is None:
            raise ValueError("EL is required for a single visit/runDir.")
        if self.runNumber is None or self.rot is None:
            raise ValueError(
                "runNumber and rot are required for a single visit/runDir."
            )

    def _fitSingleRunDir(self, runDir, arm, EL, runNumber, rot):
        return circleFit(
            runDir,
            arm,
            EL,
            runNumber=runNumber,
            rot=rot,
            xml=self.xml,
            thetaHdf5=self.thetaHdf5,
            phiHdf5=self.phiHdf5,
            update=self.update,
        )

    def _handleBatchError(self, item, exc):
        """Record a batch failure and either skip it or ask first."""

        message = f"Error while processing {item}: {exc}"

        if self.batchErrorPolicy == "prompt":
            if not askYesNo(
                f"{message}\nSkip this item and continue? [y/n] "
            ):
                raise exc

        self.errors.append(
            {
                "item": str(item),
                "error": str(exc),
                "exception": exc,
            }
        )
        warnings.warn(
            f"{message}; skipping it.",
            RuntimeWarning,
        )

    def _fitRunDirFile(self, filename):
        table = np.loadtxt(filename, dtype=str, ndmin=2)
        if table.shape[1] != 5:
            raise ValueError(
                f"{filename} must contain: "
                "runNumber, runDir, arm, EL, and rot."
            )

        rowCount = len(table)
        runNumbers = np.full(rowCount, None, dtype=object)
        runDirs = np.full(rowCount, None, dtype=object)
        arms = np.full(rowCount, None, dtype=object)
        elevations = np.full(rowCount, np.nan, dtype=float)
        rotations = np.full(rowCount, np.nan, dtype=float)

        for index, row in enumerate(table):
            try:
                runNumber = int(row[0])
                runDir = str(row[1])
                arm = str(row[2]).lower()
                if arm not in {"theta", "phi"}:
                    raise ValueError(
                        f"arm must be theta or phi, not {row[2]!r}"
                    )
                EL = normalizeEL(row[3])
                rot = normalizeEL(row[4])
                runDir2Date(runDir)
                self._fitSingleRunDir(runDir, arm, EL, runNumber, rot)
                runNumbers[index] = runNumber
                runDirs[index] = runDir
                arms[index] = arm
                elevations[index] = EL
                rotations[index] = rot
            except Exception as exc:
                self._handleBatchError(
                    f"row {index + 1} ({' '.join(row)})",
                    exc,
                )

        self.runNumbers = runNumbers
        self.runDirs = runDirs
        self.arms = arms
        self.ELs = elevations
        self.rots = rotations
        return runDirs

    def determineCenter(self):
        if self.runDir is not None:
            runDirInput = str(self.runDir)
            try:
                runDir2Date(runDirInput)
                isSingle = True
            except ValueError:
                isSingle = False

            if isSingle:
                self._validateSingleInput()
                self.runDir = self._fitSingleRunDir(
                    runDirInput,
                    self.arm,
                    self.EL,
                    self.runNumber,
                    self.rot,
                )
                return self.runDir

            return self._fitRunDirFile(runDirInput)

        if self.visit is None:
            raise ValueError("Neither visit nor runDir is defined.")

        if isinstance(self.visit, (int, np.integer)):
            self._validateSingleInput()
            self.runDir = visit2runDir(
                int(self.visit),
                self.arm,
                self.EL,
                visitRunDirList=self.visitRunDirList,
                update=self.update,
            )
            self._fitSingleRunDir(
                self.runDir,
                self.arm,
                self.EL,
                self.runNumber,
                self.rot,
            )
            return self.runDir

        if isinstance(self.visit, (str, pathlib.Path)):
            (
                self.runNumbers,
                self.runDirs,
                self.visits,
                self.arms,
                self.ELs,
                self.rots,
            ) = visitList2runDirList(
                str(self.visit),
                visitRunDirList=self.visitRunDirList,
                update=self.update,
                errorHandler=self._handleBatchError,
            )
            for runNumber, runDir, visit, arm, EL, rot in zip(
                self.runNumbers,
                self.runDirs,
                self.visits,
                self.arms,
                self.ELs,
                self.rots,
            ):
                if runDir is not None:
                    try:
                        self._fitSingleRunDir(
                            runDir,
                            arm,
                            EL,
                            runNumber,
                            rot,
                        )
                    except Exception as exc:
                        self._handleBatchError(
                            f"runDir {runDir} (arm={arm}, EL={EL})",
                            exc,
                        )
            return self.runDirs

        raise TypeError("visit must be an integer or a filename.")

    def determineLength(self):
        output = mergeGeometryByRunELRotArm(
            self.thetaHdf5,
            self.phiHdf5,
            self.geometryHdf5,
            overwrite=self.overwriteCombined,
        )
        calCombinedGeometry(self.geometryHdf5)
        return output

    def determineMovementLimits(self, outlierTolerance=1.0):
        return calculateMovementLimits(
            self.thetaHdf5,
            self.phiHdf5,
            self.geometryHdf5,
            xml=self.xml,
            outlierTolerance=outlierTolerance,
            update=self.update,
            errorHandler=self._handleBatchError,
        )

    def run(
        self,
        combine=True,
        calculateLimits=True,
        outlierTolerance=1.0,
    ):
        result = self.determineCenter()
        if combine:
            self.determineLength()
        if calculateLimits:
            self.determineMovementLimits(outlierTolerance)
        return result


geometryMeasurement = GeometryMeasurement


def parseVisit(value):
    try:
        return int(value)
    except ValueError:
        return value


def buildParser():
    parser = argparse.ArgumentParser(
        description="Fit and combine PFS cobra geometry."
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("-visit", "--visit", type=parseVisit)
    inputs.add_argument("-runDir", "--runDir")
    parser.add_argument("-arm", "--arm", choices=("theta", "phi"))
    parser.add_argument("-EL", "--EL", type=float)
    parser.add_argument("-runNumber", "--runNumber", "-run", type=int)
    parser.add_argument("-rot", "--rot", type=float)
    parser.add_argument("--xml", default=DEFAULT_XML)
    parser.add_argument(
        "--visit-runDir-list", default="visit_runDir.csv"
    )
    parser.add_argument(
        "--theta-hdf5", default="geometry/thetaGeometry.hdf5"
    )
    parser.add_argument(
        "--phi-hdf5", default="geometry/phiGeometry.hdf5"
    )
    parser.add_argument(
        "--geometry-hdf5", default="geometry/Geometry.hdf5"
    )
    parser.add_argument("--no-update", action="store_true")
    parser.add_argument("--no-combine", action="store_true")
    parser.add_argument("--no-movement-limits", action="store_true")
    parser.add_argument("--outlier-tolerance", type=float, default=1.0)
    return parser


def main():
    parser = buildParser()
    args = parser.parse_args()

    scalarInput = (
        isinstance(args.visit, int)
        or (
            args.runDir is not None
            and pathlib.Path(args.runDir).suffix == ""
        )
    )
    if scalarInput and any(
        value is None
        for value in (args.arm, args.EL, args.runNumber, args.rot)
    ):
        parser.error(
            "-arm, -EL, -runNumber, and -rot are required "
            "for a single input."
        )

    measurement = GeometryMeasurement(
        visit=args.visit,
        runDir=args.runDir,
        arm=args.arm,
        EL=args.EL,
        runNumber=args.runNumber,
        rot=args.rot,
        xml=args.xml,
        visitRunDirList=args.visit_runDir_list,
        thetaHdf5=args.theta_hdf5,
        phiHdf5=args.phi_hdf5,
        geometryHdf5=args.geometry_hdf5,
        update=not args.no_update,
        batchErrorPolicy="prompt",
    )
    measurement.run(
        combine=not args.no_combine,
        calculateLimits=not args.no_movement_limits,
        outlierTolerance=args.outlier_tolerance,
    )


if __name__ == "__main__":
    main()
