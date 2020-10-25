package ggeodesic

//Capability constants shared between Geodesic and GeodesicLine.

const capNone = 0
const capC1 = 1 << 0
const capC1p = 1 << 1
const capC2 = 1 << 2
const capC3 = 1 << 3
const capC4 = 1 << 4
const capAll = 0x1F
const capMask = capAll
const outAll = 0x7F80
const outMask = 0xFF80 // Includes LONG_UNROLL
const empty = 0
const LATITUDE = 1<<7 | capNone
const LONGITUDE = 1<<8 | capC3
const AZIMUTH = 1<<9 | capNone
const DISTANCE = 1<<10 | capC1
const STANDARD = LATITUDE | LONGITUDE | AZIMUTH | DISTANCE
const distanceIn = 1<<11 | capC1 | capC1p
const reducedLength = 1<<12 | capC1 | capC2
const geodesicScale = 1<<13 | capC1 | capC2
const area = 1<<14 | capC4
const longUnroll = 1 << 15
const all = outAll | capAll // Does not include LONG_UNROLL
