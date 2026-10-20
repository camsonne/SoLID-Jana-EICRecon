#!/usr/bin/env python3
"""Convert the SoLID CLEO-II magnet of solid_gemc (GEMC text geometry) to a DD4hep compact file.

Inputs (from https://github.com/JeffersonLab/solid_gemc, geometry/magnet_moved, the "moved" /
longer-endcap layout used with the DDVCS GEM positions):
    solid_solenoid__geometry_v2.txt        CLEO iron: yoke octagons, coil collars, upstream plug,
                                           endcap nose cone, endcap cylinder and plates, corrector
    solid_solenoid_cryo__geometry_Original.txt   cryostat and coil
Every GEMC volume is placed at the origin, as in solid_gemc. GEMC Pgon (octagonal G4Polyhedra,
radii to the flat faces) volumes become DD4hep_PolyconeSupport polycones with the same radii, i.e.
about 5% less iron than the real octagon but no overlap with detectors near its corners; Polycone
and Tube volumes are converted exactly.

usage: gemc_magnet_to_dd4hep.py yoke.txt cryo.txt > geometry/compact/magnet_ddvcs.xml
"""
import sys

MATERIALS = {"G4_Fe": "Iron", "G4_STAINLESS-STEEL": "StainlessSteel", "G4_Al": "Aluminum", "G4_Cu": "Copper"}


def num(tok):
    tok = tok.strip()
    for unit in ("*cm", "*deg", "*counts", "*mm"):
        if tok.endswith(unit):
            v = float(tok[: -len(unit)])
            return v / 10.0 if unit == "*mm" else v
    return float(tok)


def volumes(path):
    for line in open(path):
        f = [x.strip() for x in line.split("|")]
        if len(f) < 9 or not f[0]:
            continue
        name, pos, shape, dims, mat = f[0], f[3], f[6], f[7].split(), f[8]
        z0 = num(pos.split()[2])
        if shape in ("Pgon", "Polycone"):
            if shape == "Pgon":
                n = int(num(dims[3]))
                vals = [num(x) for x in dims[4:]]
            else:
                n = int(num(dims[2]))
                vals = [num(x) for x in dims[3:]]
            rin, rout, z = vals[:n], vals[n:2 * n], vals[2 * n:3 * n]
        elif shape == "Tube":
            rmin, rmax, hz = (num(x) for x in dims[:3])
            rin, rout, z = [rmin, rmin], [rmax, rmax], [-hz, hz]
        else:
            raise ValueError(f"{name}: unsupported shape {shape}")
        yield name, shape, rin, rout, [zz + z0 for zz in z], MATERIALS[mat]


def main():
    out = ['<lccdd>',
           '  <!-- SoLID CLEO-II solenoid: coil, cryostat and iron (solid_gemc geometry/magnet_moved,',
           '       solid_solenoid__geometry_v2 + solid_solenoid_cryo), converted by scripts/gemc_magnet_to_dd4hep.py.',
           '       Octagonal pieces are polycones at the octagon flat-face radius (~5% less iron). -->',
           '  <display>',
           '    <vis name="MagnetIronVis" r="0.6" g="0.25" b="0.95" alpha="0.3" showDaughters="false" visible="true"/>',
           '    <vis name="MagnetCoilVis" r="1.0" g="0.5" b="0.07" alpha="0.5" showDaughters="false" visible="true"/>',
           '  </display>',
           '  <detectors>']
    det_id = 80
    for path in sys.argv[1:]:
        for name, shape, rin, rout, z, mat in volumes(path):
            short = name.replace("solid_solenoid_cryo_", "Magnet").replace("solid_solenoid_", "Magnet")
            vis = "MagnetCoilVis" if mat in ("Aluminum", "Copper") else "MagnetIronVis"
            out.append(f'    <!-- {name} ({shape}) -->')
            out.append(f'    <detector id="{det_id}" name="{short}" type="DD4hep_PolyconeSupport" vis="{vis}">')
            out.append(f'      <material name="{mat}"/>')
            for a, b, c in zip(rin, rout, z):
                out.append(f'      <zplane rmin="{a:g}*cm" rmax="{b:g}*cm" z="{c:g}*cm"/>')
            out.append('    </detector>')
            det_id += 1
    out += ['  </detectors>', '</lccdd>']
    print("\n".join(out))


if __name__ == "__main__":
    main()
