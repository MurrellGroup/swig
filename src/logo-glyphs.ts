/**
 * Literal glyph contours from DejaVu Sans Mono Bold.
 *
 * Coordinates retain the TrueType convention (positive y points upward).
 * Copyright (c) 2003 by Bitstream, Inc. DejaVu changes are public domain.
 * Distributed under the Bitstream Vera license in DEJAVU_FONT_LICENSE.txt.
 */
export const LOGO_GLYPH_FONT = "DejaVu Sans Mono Bold";

export interface LogoGlyphOutline {
  readonly d: string;
  readonly xMin: number;
  readonly yMin: number;
  readonly xMax: number;
  readonly yMax: number;
  readonly advance: number;
}

export const LOGO_GLYPH_OUTLINES = {
  "A": { d: "M616 1223 477 612H756ZM436 1493H797L1200 0H905L813 369H418L328 0H33Z", xMin: 33, yMin: 0, xMax: 1200, yMax: 1493, advance: 1233 },
  "C": { d: "M1081 43Q1011 7 934.0 -11.0Q857 -29 772 -29Q470 -29 311.0 170.0Q152 369 152 745Q152 1122 311.0 1321.0Q470 1520 772 1520Q857 1520 935.0 1502.0Q1013 1484 1081 1448V1120Q1005 1190 933.5 1222.5Q862 1255 786 1255Q624 1255 541.5 1126.5Q459 998 459 745Q459 493 541.5 364.5Q624 236 786 236Q862 236 933.5 268.5Q1005 301 1081 371Z", xMin: 152, yMin: -29, xMax: 1081, yMax: 1520, advance: 1233 },
  "D": { d: "M432 1227V266H512Q686 266 760.0 375.5Q834 485 834 748Q834 1009 760.0 1118.0Q686 1227 512 1227ZM137 1493H453Q819 1493 980.0 1318.5Q1141 1144 1141 748Q1141 351 980.0 175.5Q819 0 453 0H137Z", xMin: 137, yMin: 0, xMax: 1141, yMax: 1493, advance: 1233 },
  "E": { d: "M1098 0H168V1493H1098V1233H463V911H1038V651H463V260H1098Z", xMin: 168, yMin: 0, xMax: 1098, yMax: 1493, advance: 1233 },
  "F": { d: "M1112 1233H477V911H1055V651H477V0H182V1493H1112Z", xMin: 182, yMin: 0, xMax: 1112, yMax: 1493, advance: 1233 },
  "G": { d: "M872 270V555H670V803H1130V119Q1045 46 942.5 8.5Q840 -29 723 -29Q433 -29 275.0 172.5Q117 374 117 745Q117 1122 276.5 1321.0Q436 1520 737 1520Q827 1520 914.0 1494.5Q1001 1469 1077 1421V1094Q1015 1174 934.5 1214.5Q854 1255 758 1255Q590 1255 507.0 1128.5Q424 1002 424 745Q424 496 504.0 366.0Q584 236 737 236Q783 236 817.0 244.5Q851 253 872 270Z", xMin: 117, yMin: -29, xMax: 1130, yMax: 1520, advance: 1233 },
  "H": { d: "M137 1493H432V924H801V1493H1096V0H801V664H432V0H137Z", xMin: 137, yMin: 0, xMax: 1096, yMax: 1493, advance: 1233 },
  "I": { d: "M172 1233V1493H1061V1233H764V260H1061V0H172V260H469V1233Z", xMin: 172, yMin: 0, xMax: 1061, yMax: 1493, advance: 1233 },
  "K": { d: "M117 1493H412V903L874 1493H1208L737 905L1225 0H897L543 672L412 506V0H117Z", xMin: 117, yMin: 0, xMax: 1225, yMax: 1493, advance: 1233 },
  "L": { d: "M225 0V1493H520V260H1151V0Z", xMin: 225, yMin: 0, xMax: 1151, yMax: 1493, advance: 1233 },
  "M": { d: "M86 1493H438L616 838L793 1493H1147V0H893V1196L735 543H500L340 1196V0H86Z", xMin: 86, yMin: 0, xMax: 1147, yMax: 1493, advance: 1233 },
  "N": { d: "M119 1493H436L852 408V1493H1112V0H797L379 1085V0H119Z", xMin: 119, yMin: 0, xMax: 1112, yMax: 1493, advance: 1233 },
  "P": { d: "M457 1245V807H578Q723 807 781.5 856.0Q840 905 840 1026Q840 1147 781.5 1196.0Q723 1245 578 1245ZM162 1493H567Q876 1493 1011.5 1383.0Q1147 1273 1147 1026Q1147 779 1011.5 669.0Q876 559 567 559H457V0H162Z", xMin: 162, yMin: 0, xMax: 1147, yMax: 1493, advance: 1233 },
  "Q": { d: "M656 -23Q642 -26 632.5 -27.5Q623 -29 614 -29Q357 -29 224.5 167.0Q92 363 92 745Q92 1128 224.5 1324.0Q357 1520 616 1520Q876 1520 1008.5 1324.0Q1141 1128 1141 745Q1141 482 1078.0 304.5Q1015 127 895 51L1081 -131L879 -281ZM616 1255Q503 1255 451.0 1134.5Q399 1014 399 745Q399 477 451.0 356.5Q503 236 616 236Q730 236 782.0 356.5Q834 477 834 745Q834 1014 782.0 1134.5Q730 1255 616 1255Z", xMin: 92, yMin: -281, xMax: 1141, yMax: 1520, advance: 1233 },
  "R": { d: "M807 705Q851 696 883.5 663.5Q916 631 963 537L1233 0H909L729 377Q721 393 708 421Q629 590 522 590H428V0H133V1493H559Q847 1493 972.5 1391.0Q1098 1289 1098 1059Q1098 905 1023.0 814.0Q948 723 807 705ZM428 1245V838H567Q688 838 740.5 885.5Q793 933 793 1042Q793 1151 741.0 1198.0Q689 1245 567 1245Z", xMin: 133, yMin: 0, xMax: 1233, yMax: 1493, advance: 1233 },
  "S": { d: "M510 655Q287 740 208.0 833.5Q129 927 129 1085Q129 1288 259.0 1404.0Q389 1520 616 1520Q719 1520 822.0 1496.5Q925 1473 1026 1427V1139Q931 1206 833.0 1241.0Q735 1276 639 1276Q532 1276 475.0 1233.0Q418 1190 418 1110Q418 1048 459.5 1007.5Q501 967 633 918L760 870Q940 804 1025.0 695.0Q1110 586 1110 420Q1110 194 976.5 82.5Q843 -29 573 -29Q462 -29 350.5 -2.5Q239 24 135 76V381Q253 297 363.5 256.0Q474 215 582 215Q691 215 751.0 264.5Q811 314 811 403Q811 470 771.0 520.5Q731 571 655 600Z", xMin: 129, yMin: -29, xMax: 1110, yMax: 1520, advance: 1233 },
  "T": { d: "M764 0H469V1235H90V1493H1143V1235H764Z", xMin: 90, yMin: 0, xMax: 1143, yMax: 1493, advance: 1233 },
  "V": { d: "M616 246 879 1493H1176L821 0H412L57 1493H354Z", xMin: 57, yMin: 0, xMax: 1176, yMax: 1493, advance: 1233 },
  "W": { d: "M0 1493H258L365 397L494 1106H739L889 397L973 1493H1233L1061 0H786L616 784L457 0H184Z", xMin: 0, yMin: 0, xMax: 1233, yMax: 1493, advance: 1233 },
  "Y": { d: "M8 1493H326L616 893L907 1493H1225L764 588V0H469V588Z", xMin: 8, yMin: 0, xMax: 1225, yMax: 1493, advance: 1233 },
  "*": { d: "M1108 1217 778 1044 1108 870 1032 729 700 913V569H528V913L197 729L121 870L453 1044L121 1217L197 1358L528 1176V1520H700V1176L1032 1358Z", xMin: 121, yMin: 569, xMax: 1108, yMax: 1520, advance: 1233 },
  "X": { d: "M1206 0H901L616 494L332 0H27L465 758L39 1493H344L616 1018L889 1493H1194L770 758Z", xMin: 27, yMin: 0, xMax: 1206, yMax: 1493, advance: 1233 },
  "-": { d: "M301 735H932V444H301Z", xMin: 301, yMin: 444, xMax: 932, yMax: 735, advance: 1233 },
} as const satisfies Readonly<Record<string, LogoGlyphOutline>>;

export interface LogoGlyphPathPlacement {
  readonly d: string;
  readonly x: number;
}

export interface LogoGlyphRun {
  readonly paths: readonly LogoGlyphPathPlacement[];
  readonly xMin: number;
  readonly yMin: number;
  readonly xMax: number;
  readonly yMax: number;
  readonly width: number;
  readonly height: number;
}

const runCache = new Map<string, LogoGlyphRun | undefined>();

/** Compose one or more literal outlines into a fixed-advance glyph run. */
export function logoGlyphRun(symbol: string): LogoGlyphRun | undefined {
  if (runCache.has(symbol)) return runCache.get(symbol);
  let cursor = 0;
  let xMin = Number.POSITIVE_INFINITY;
  let yMin = Number.POSITIVE_INFINITY;
  let xMax = Number.NEGATIVE_INFINITY;
  let yMax = Number.NEGATIVE_INFINITY;
  const paths: LogoGlyphPathPlacement[] = [];
  for (const character of symbol) {
    const outline = LOGO_GLYPH_OUTLINES[character as keyof typeof LOGO_GLYPH_OUTLINES];
    if (!outline) {
      runCache.set(symbol, undefined);
      return undefined;
    }
    paths.push({ d: outline.d, x: cursor });
    xMin = Math.min(xMin, cursor + outline.xMin);
    yMin = Math.min(yMin, outline.yMin);
    xMax = Math.max(xMax, cursor + outline.xMax);
    yMax = Math.max(yMax, outline.yMax);
    cursor += outline.advance;
  }
  if (!paths.length || !(xMax > xMin) || !(yMax > yMin)) {
    runCache.set(symbol, undefined);
    return undefined;
  }
  const run: LogoGlyphRun = { paths, xMin, yMin, xMax, yMax, width: xMax - xMin, height: yMax - yMin };
  runCache.set(symbol, run);
  return run;
}
