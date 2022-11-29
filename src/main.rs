use std::{
    env,
    fs::File,
    io::{BufRead, BufReader, Write},
};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        panic!("Parameters count error!");
    }
    let test: [char; 7] = ['A', 'T', 'O', 'M', ' ', ' ', '\0'];
    let mut line2: [char; 7] = [' '; 7];
    const SIZE: usize = 250000;
    const BXMX: usize = 200000;
    let mut atmnum: usize;
    let mut name_temp: char;
    let mut name_temp2: [char; 4] = [' '; 4];
    let mut name: Vec<usize> = Vec::with_capacity(SIZE);
    name.resize(SIZE, 0);
    let mut bnam: Vec<i32> = Vec::with_capacity(SIZE);
    bnam.resize(SIZE, 0);
    let mut alt_loc: char;
    let mut res_name: [char; 4] = [' '; 4];
    let mut chain_id: Vec<char> = Vec::with_capacity(SIZE);
    chain_id.resize(SIZE, ' ');
    let mut res_seq_temp: [char; 5] = [' '; 5];
    let mut res_seq: Vec<i32> = Vec::with_capacity(SIZE);
    res_seq.resize(SIZE, 0);
    let mut resnum: Vec<i32> = Vec::with_capacity(SIZE);
    resnum.resize(SIZE, 0);
    let mut x: [char; 9] = [' '; 9];
    let mut y: [char; 9] = [' '; 9];
    let mut z: [char; 9] = [' '; 9];
    let mut xyz: Vec<[f64; 3]> = Vec::with_capacity(SIZE);
    xyz.resize(SIZE, [0.0; 3]);
    let mut flag = 0;
    let mut kadd: i32;
    const CHAINDIF: i32 = 10000;
    let mut min: [f64; 4] = [0.0; 4];
    let mut max: [f64; 4] = [0.0; 4];
    let mut nbx: [i32; 4] = [0; 4];
    const BOXSIZE: f64 = 4.0;
    let mut ibox1: Vec<[i32; 16]> = Vec::with_capacity(BXMX);
    ibox1.resize(BXMX, [0; 16]);
    let mut temp: usize;
    let mut ix: i32;
    let mut iy: i32;
    let mut iz: i32;
    let mut ind: usize;
    let mut most = 0;
    const RADIUS: f64 = 3.75;
    const RSQ: f64 = RADIUS * RADIUS;
    const RADMIN: f64 = 3.25;
    const SSQ: f64 = RADMIN * RADMIN;
    let ndelta: i32 = (RADIUS / BOXSIZE).ceil() as i32;
    let mut dsq: f64;
    let mut jbx: i32;
    let mut jby: i32;
    let mut jbz: i32;
    let mut ibz1: i32;
    let mut ibz2: i32;
    let mut iby1: i32;
    let mut iby2: i32;
    let mut ibx1: i32;
    let mut ibx2: i32;
    let mut temp1: f64;
    let mut temp2: f64;
    let mut count = 0;
    const MAXWIN: f64 = 100.694;
    let mut c: [[f64; 4]; 4] = [[0.0; 4]; 4];
    let mut matrix: [f64; 6] = [0.0; 6];
    let mut mtrx: f64;
    let mut mtrxstat: f64 = 0.0;
    let mut pstat: f64 = 0.0;
    let mut stat: f64 = 0.0;
    let mut errat: Vec<f64> = Vec::with_capacity(SIZE);
    errat.resize(SIZE, 0.0);
    let chainx: i32;
    let mut ir1: [i32; 100] = [0; 100];
    let mut ir2: [i32; 100] = [0; 100];
    let mut id_by_chain: [char; 100] = [' '; 100];
    let mut ms: f64;
    let mut mst: f64;
    let mut np: i32;
    let mut ir0: f64;
    let mut ir: f64;
    let mut z2: usize;

    let filepath = args[1].clone();
    let logfilename = filepath.clone();
    let filepath = filepath + ".pdb";
    let logfilename = logfilename + ".logf";

    let mut fout = File::create(logfilename).expect("Failed opening errat.logf");

    let mut lmt: [f64; 3] = [0.0; 3];
    lmt[2] = 11.526684477428809;
    lmt[1] = 17.190823041860433;
    atmnum = 0;
    resnum[0] = 0;

    let mut flag2 = 0;
    kadd = 0;

    let file = File::open(&filepath).expect("Failed opening file.pdb");
    let fin = BufReader::new(file);
    let mut i = 0;
    for line in fin.lines() {
        let line = line.unwrap().chars().collect::<Vec<char>>();
        if line[10] == 'D' && line[11] == 'N' && line[12] == 'A' {
            panic!("This is a DNA")
        }
        for j in 0..6 {
            line2[j] = line[j];
        }
        line2[6] = '\0';
        if line2 == test {
            i += 1;
            if i > SIZE - 1 {
                fout.write_all(
                    "ERROR: PDB WITH TOO MANY ATOMS. CUT OFF FURTHER INPUT.\n".as_bytes(),
                )
                .expect("write failed");
                break;
            }
            name_temp = line[13];
            match name_temp {
                'C' => name[i] = 1,
                'N' => name[i] = 2,
                'O' => name[i] = 3,
                _ => name[i] = 0,
            }
            name_temp2[0] = line[13];
            name_temp2[1] = line[14];
            name_temp2[2] = line[15];
            name_temp2[3] = '\0';
            if name_temp2 == ['N', ' ', ' ', '\0'] || name_temp2 == ['C', ' ', ' ', '\0'] {
                bnam[i] = 1;
            } else {
                bnam[i] = 0;
            }
            alt_loc = line[16];
            for j in 17..20 {
                res_name[j - 17] = line[j];
            }
            res_name[3] = '\0';
            chain_id[i] = line[21];
            for j in 22..26 {
                res_seq_temp[j - 22] = line[j];
            }
            res_seq[i] = res_seq_temp
                .iter()
                .collect::<String>()
                .trim()
                .parse::<i32>()
                .unwrap();
            for j in 30..38 {
                x[j - 30] = line[j];
            }
            xyz[i][0] = x.iter().collect::<String>().trim().parse::<f64>().unwrap();
            for j in 38..46 {
                y[j - 38] = line[j];
            }
            xyz[i][1] = y.iter().collect::<String>().trim().parse::<f64>().unwrap();
            for j in 46..54 {
                z[j - 46] = line[j];
            }
            xyz[i][2] = z.iter().collect::<String>().trim().parse::<f64>().unwrap();
            if !((alt_loc == ' ') || (alt_loc == 'A') || (alt_loc == 'a') || (alt_loc == 'P')) {
                fout.write_all(
                    format!("Reject 2' Conformation atom#	{}	chain	{}\n", i, chain_id[i]).as_bytes(),
                )
                .expect("write failed");
                i -= 1;
                flag = 1;
            }
            if !(res_name == ['G', 'L', 'Y', '\0']
                || res_name == ['A', 'L', 'A', '\0']
                || res_name == ['V', 'A', 'L', '\0']
                || res_name == ['L', 'E', 'U', '\0']
                || res_name == ['I', 'L', 'E', '\0']
                || res_name == ['T', 'Y', 'R', '\0']
                || res_name == ['C', 'Y', 'S', '\0']
                || res_name == ['M', 'E', 'T', '\0']
                || res_name == ['T', 'R', 'P', '\0']
                || res_name == ['P', 'H', 'E', '\0']
                || res_name == ['H', 'I', 'S', '\0']
                || res_name == ['P', 'R', 'O', '\0']
                || res_name == ['S', 'E', 'R', '\0']
                || res_name == ['T', 'H', 'R', '\0']
                || res_name == ['L', 'Y', 'S', '\0']
                || res_name == ['A', 'R', 'G', '\0']
                || res_name == ['G', 'L', 'U', '\0']
                || res_name == ['A', 'S', 'P', '\0']
                || res_name == ['G', 'L', 'N', '\0']
                || res_name == ['A', 'S', 'N', '\0'])
            {
                i -= 1;
                flag = 1;
                fout.write_all(
                    format!("***Warning: Reject Nonstardard Residue - {:?}\n", res_name).as_bytes(),
                )
                .expect("write failed");
            }
            if (chain_id[i] != chain_id[i - 1]) && (i >= 2) && (flag != 1) {
                kadd += 1;
                fout.write_all(format!("INCREMENTING CHAIN (kadd) {}\n", kadd).as_bytes())
                    .expect("write failed");
            }
            if flag != 1 {
                resnum[i] = res_seq[i] + (kadd * CHAINDIF);
                atmnum = i;
            }
            if resnum[i] < resnum[i - 1] && chain_id[i] == chain_id[i - 1] && flag == 0 && i >= 2 {
                fout.write_all(
                    format!(
                        "ERROR: RESNUM DECREASE. TERMINATE ANALYSIS{}  {}\n",
                        resnum[i],
                        resnum[i - 1]
                    )
                    .as_bytes(),
                )
                .expect("write failed");
            }
            if i > 2
                && resnum[i] != resnum[i - 1]
                && chain_id[i] == chain_id[i - 1]
                && flag == 0
                && (resnum[i] - resnum[i - 1]) > 1
            {
                fout.write_all(
                    format!(
                        "WARNING: Missing Residues{}>>>{}\n",
                        resnum[i - 1],
                        resnum[i]
                    )
                    .as_bytes(),
                )
                .expect("write failed");
            }
            errat[(resnum[i] + 4) as usize] = 0.0;
            flag = 0;
        }
    }
    for i in 1..=3 {
        min[i] = 999.0;
        max[i] = -999.0;
    }
    for i in 1..=atmnum {
        for j in 1..=3 {
            if xyz[i][j - 1] < min[j] {
                min[j] = xyz[i][j - 1];
            }
            if xyz[i][j - 1] > max[j] {
                max[j] = xyz[i][j - 1];
            }
        }
    }
    for j in 1..=3 {
        fout.write_all(format!("{} ", min[j]).as_bytes())
            .expect("write failed");
    }
    for j in 1..=3 {
        fout.write_all(format!("{} ", max[j]).as_bytes())
            .expect("write failed");
    }
    fout.write_all("\n".as_bytes()).expect("write failed");
    for i in 1..=3 {
        nbx[i] = ((max[i] - min[i]) / BOXSIZE) as i32 + 1;
    }
    if nbx[1] * nbx[2] * nbx[3] > BXMX as i32 - 1 {
        fout.write_all("ERROR: TOO MANY BOXES\n".as_bytes())
            .expect("write failed");
        flag2 = 1;
    }
    if flag2 != 1 {
        for i in 1..=atmnum {
            ix = ((xyz[i][0] - (min[1] - 0.00001)) / BOXSIZE) as i32;
            iy = ((xyz[i][1] - (min[2] - 0.00001)) / BOXSIZE) as i32;
            iz = ((xyz[i][2] - (min[3] - 0.00001)) / BOXSIZE) as i32;
            ind = (1 + ix + iy * nbx[1] + iz * nbx[1] * nbx[2]) as usize;

            ibox1[ind][0] = ibox1[ind][0] + 1;
            temp = ibox1[ind][0] as usize;
            if temp < 16 {
                ibox1[ind][temp] = i as i32;
            }
        }
        for i in 1..=nbx[1] * nbx[2] * nbx[3] {
            let i = i as usize;
            if ibox1[i][0] > 15 {
                fout.write_all(format!("TOO MANY ATOMS IN BOX #:	{}\n", ibox1[0][i]).as_bytes())
                    .expect("write failed");
                let mut j = 1;
                while j < ibox1[i][0] && ibox1[i][0] < 16 {
                    fout.write_all(format!("ibox[{}][{}]", j, i).as_bytes())
                        .expect("write failed");
                    j += 1;
                }
                flag2 = 1;
            }
            if ibox1[i][0] > most {
                most = ibox1[i][0];
            }
        }
        pstat = 0.0;
        stat = 0.0;
        mtrxstat = 0.0;
        if flag2 != 1 {
            for i in 1..=atmnum {
                if resnum[i] > resnum[i - 1] || i == 1 {
                    for aa in 0..4 {
                        for ab in 0..4 {
                            c[ab][aa] = 0.0;
                        }
                    }
                    let mut s = 1;
                    let mut v = i;
                    while s < 10 && v <= atmnum {
                        if (resnum[v + 1] - resnum[v] < 100 && (resnum[v + 1] - resnum[v] > 0))
                            || v == atmnum
                        {
                            s += 1;
                        }
                        v += 1;
                    }
                    v -= 1;
                    if res_seq[v] > res_seq[i] && s == 10 {
                        for rer in i..=v {
                            jbx = ((xyz[rer][0] - (min[1] - 0.00001)) / BOXSIZE) as i32;
                            jby = ((xyz[rer][1] - (min[2] - 0.00001)) / BOXSIZE) as i32;
                            jbz = ((xyz[rer][2] - (min[3] - 0.00001)) / BOXSIZE) as i32;
                            ibz1 = jbz - ndelta;
                            if ibz1 < 0 {
                                ibz1 = 0;
                            }
                            ibz2 = jbz + ndelta;
                            if ibz2 > nbx[3] - 1 {
                                ibz2 = nbx[3] - 1;
                            }
                            iby1 = jby - ndelta;
                            if iby1 < 0 {
                                iby1 = 0;
                            }
                            iby2 = jby + ndelta;
                            if iby2 > nbx[2] - 1 {
                                iby2 = nbx[2] - 1;
                            }
                            ibx1 = jbx - ndelta;
                            if ibx1 < 0 {
                                ibx1 = 0;
                            }
                            ibx2 = jbx + ndelta;
                            if ibx2 > nbx[1] - 1 {
                                ibx2 = nbx[1] - 1;
                            }
                            for j in ibz1..=ibz2 {
                                for k in iby1..=iby2 {
                                    for l in ibx1..=ibx2 {
                                        ind = (1 + l + k * nbx[1] + j * nbx[1] * nbx[2]) as usize;
                                        for m in 1..=ibox1[ind][0] {
                                            let n: usize = ibox1[ind][m as usize] as usize;
                                            if resnum[rer] != resnum[n] {
                                                dsq = 0.0;
                                                for p in 0..=2 {
                                                    dsq = dsq + (xyz[n][p] - xyz[rer][p]).powf(2.0);
                                                }
                                                if dsq < RSQ {
                                                    if !(bnam[rer] == 1
                                                        && bnam[n] == 1
                                                        && (((resnum[n] == resnum[rer] + 1)
                                                            && (name[rer] == 1)
                                                            && (name[n] == 2))
                                                            || ((resnum[rer] == resnum[n] + 1)
                                                                && (name[rer] == 2)
                                                                && (name[n] == 1))))
                                                    {
                                                        if n >= i && n <= v {
                                                            if resnum[rer] > resnum[n] {
                                                                if dsq <= SSQ {
                                                                    temp1 = 1.0;
                                                                } else {
                                                                    temp1 =
                                                                        2.0 * (3.75 - dsq.sqrt());
                                                                }
                                                                count += 1;
                                                                c[name[n]][name[rer]] =
                                                                    c[name[n]][name[rer]] + temp1;
                                                            }
                                                        } else {
                                                            if dsq <= SSQ {
                                                                temp1 = 1.0;
                                                            } else {
                                                                temp1 = 2.0 * (3.75 - dsq.sqrt());
                                                            }
                                                            count += 1;
                                                            c[name[n]][name[rer]] =
                                                                c[name[n]][name[rer]] + temp1;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        temp2 = 0.0;
                        for q in 1..=3 {
                            for r in 1..=3 {
                                temp2 += c[r][q];
                            }
                        }
                        if temp2 <= 0.0 {
                            fout.write_all(
                                format!(
                                    "{} {} {} WARNING: No Interactions in This Frame\n",
                                    temp2,
                                    resnum[i] + 4,
                                    count
                                )
                                .as_bytes(),
                            )
                            .expect("write failed");
                        }
                        if temp2 > MAXWIN {
                            matrix[1] = c[1][1] / temp2;
                            matrix[2] = (c[2][1] + c[1][2]) / temp2;
                            matrix[3] = (c[3][1] + c[1][3]) / temp2;
                            matrix[4] = c[2][2] / temp2;
                            matrix[5] = (c[3][2] + c[2][3]) / temp2;
                            mtrx = matrixdb(matrix);
                            stat += 1.0;
                            mtrxstat += mtrx;
                            if mtrx > lmt[1] {
                                pstat += 1.0;
                            } else if mtrx > lmt[2] {
                                pstat += 1.0;
                            }
                            errat[(resnum[i] + 4) as usize] = mtrx;
                        } else {
                            fout.write_all(
                                format!(
                                    "WARNING: Frame	{} Below Minimum Interaction Limit.\n",
                                    resnum[i] + 4
                                )
                                .as_bytes(),
                            )
                            .expect("write failed");
                        }
                    }
                }
            }
        }
    }

    if stat > 0.0 {
        fout.write_all(
            format!(
                "Total frames: {} P frames: {} Number: {}.\n",
                stat,
                pstat,
                pstat / stat
            )
            .as_bytes(),
        )
        .expect("write failed");
        fout.write_all(format!("Avg Probability {}\n", mtrxstat / stat).as_bytes())
            .expect("write failed");
        fout.write_all(
            format!(
                "# Overall quality factor: {}\n",
                100.0 - (100.0 * pstat / stat)
            )
            .as_bytes(),
        )
        .expect("write failed");
        print!("{}", 100.0 - (100.0 * pstat / stat));
        chainx = 1 + (resnum[atmnum] - 4) / 10000;
        z2 = 1;
        ir1[z2] = resnum[1] + 4;
        ir2[z2] = 0;
        id_by_chain[z2] = chain_id[1];
        for z1 in 1..=atmnum {
            if z1 == atmnum - 1 {
                ir2[z2] = resnum[atmnum] - 4;
            } else if chain_id[z1] != chain_id[z1 + 1] && resnum[z1] > 4 {
                ir2[z2] = resnum[z1] - 4;
                z2 += 1;
                ir1[z2] = resnum[z1 + 1] + 4;
                id_by_chain[z2] = chain_id[z1 + 1];
            }
        }
        mst = 0.0;
        for ich in 1..=chainx {
            let ich = ich as usize;
            ms = (ir2[ich] - ir1[ich] + 1) as f64 / 301.0;
            ms = (ir2[ich] - ir1[ich] + 1) as f64 / ms;
            if ms > mst {
                mst = ms;
            }
            if mst < 200.0 {
                mst = 200.0;
            }
        }
        for ich in 1..=chainx {
            let ich = ich as usize;
            np = 1 + ((ir2[ich] - ir1[ich] + 1) / mst as i32);
            for z1 in 1..=np {
                ir0 = ir1[ich] as f64 + mst * (z1 - 1) as f64;
                ir = ir0 + mst - 1.0;
                if ir > ir2[ich] as f64 {
                    ir = ir2[ich] as f64;
                }
                fout.write_all(
                    format!(
                        "# Chain Label  {}: Residue range {} to {}\n",
                        id_by_chain[ich], ir0, ir
                    )
                    .as_bytes(),
                )
                .expect("write failed");
            }
        }
    }
}

fn matrixdb(mut matrix: [f64; 6]) -> f64 {
    let mut c1: [[f64; 6]; 6] = [[0.0; 6]; 6];
    let mut d1: [[f64; 6]; 6] = [[0.0; 6]; 6];
    let avg: [f64; 6] = [
        0.0,
        0.192765509919262,
        0.195575208778518,
        0.275322406824210,
        0.059102357035642,
        0.233154192767480,
    ];
    let b1: [[f64; 6]; 6] = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [
            0.0,
            5040.279078850848200,
            3408.805141583649400,
            4152.904423767300600,
            4236.200004171890200,
            5054.781210204625500,
        ],
        [
            0.0,
            3408.805141583648900,
            8491.906094010220800,
            5958.881777877950300,
            1521.387352718486200,
            4304.078200827221700,
        ],
        [
            0.0,
            4152.904423767301500,
            5958.881777877952100,
            7637.167089335050100,
            6620.715738223072500,
            5287.691183798410700,
        ],
        [
            0.0,
            4236.200004171890200,
            1521.387352718486200,
            6620.715738223072500,
            18368.343774298410000,
            4050.797811118806700,
        ],
        [
            0.0,
            5054.781210204625500,
            4304.078200827220800,
            5287.691183798409800,
            4050.797811118806700,
            6666.856740479164700,
        ],
    ];
    let mut a: [[f64; 6]; 6] = [[0.0; 6]; 6];
    let mut b: [[f64; 6]; 6] = [[0.0; 6]; 6];
    // fout.write_all("PROCESSING FRAME STATISTICS:\n".as_bytes())
    //     .expect("write failed");

    for u in 1..6 {
        matrix[u] = matrix[u] - avg[u];
    }

    // fout.write_all("Vertical Matrix A".as_bytes())
    //     .expect("write failed");
    let m1 = 1;
    let p1 = 5;
    for u in 1..=m1 {
        // fout.write_all("\n".as_bytes()).expect("write failed");
        for v in 1..=p1 {
            a[v][u] = matrix[v];
            b[u][v] = matrix[v];
            // fout.write_all(format!("{}      ", a[v][u]).as_bytes())
            //     .expect("write failed");
        }
    }
    let n1 = 5;
    // fout.write_all(format!("\nMatrix Product AxB1 = C1 {} x {}", m1, n1).as_bytes())
    //     .expect("write failed");
    for u in 1..=m1 {
        for v in 1..=n1 {
            let mut x: f64 = 0.0;
            for k1 in 1..=p1 {
                x = x + a[k1][u] * b1[v][k1];
            }
            c1[v][u] = x;
            // fout.write_all(format!("{}      ", c1[v][u]).as_bytes())
            //     .expect("write failed");
        }
    }

    // fout.write_all("\nHorizontal Matrix B".as_bytes())
    //     .expect("write failed");
    let p1 = 5;
    let n1 = 1;
    // for u in 1..=p1 {
    //     fout.write_all("\n".as_bytes()).expect("write failed");
    //     for v in 1..=n1 {
    //         fout.write_all(format!("{}      ", b[v][u]).as_bytes())
    //             .expect("write failed");
    //     }
    // }
    for u in 1..=m1 {
        for v in 1..=n1 {
            let mut x: f64 = 0.0;
            for k1 in 1..=p1 {
                x = x + c1[k1][u] * b[v][k1];
            }
            d1[v][u] = x;
        }
    }
    // fout.write_all(format!("\nTotal Matrix\n{}      \n", d1[1][1]).as_bytes())
    //     .expect("write failed");
    return d1[1][1];
}
