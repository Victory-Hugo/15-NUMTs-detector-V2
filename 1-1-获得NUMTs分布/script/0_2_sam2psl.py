#!/usr/bin/env python3
import sys, re

if len(sys.argv) != 4:
    print("Usage: sam2psl_manual.py input.sam ref.fai output.psl")
    sys.exit(1)

sam_path, fai_path, psl_path = sys.argv[1:]

# 1) 读染色体长度
chrom_size = {}
with open(fai_path) as f:
    for line in f:
        chrom, length, *_ = line.strip().split('\t')
        chrom_size[chrom] = int(length)

# 2) CIGAR 解析正则
cigar_pat = re.compile(r'(\d+)([MIDNSHP=X])')

def parse_cigar(cigar, pos1, nm):
    q_ofs = 0
    t_ofs = pos1 - 1
    blocks=[]; qS=[]; tS=[]
    qGaps=qGapBp=tGaps=tGapBp=0
    for length,op in cigar_pat.findall(cigar):
        L=int(length)
        if op in ('M','X','='):
            blocks.append(L); qS.append(q_ofs); tS.append(t_ofs)
            q_ofs += L; t_ofs += L
        elif op=='I':
            qGaps+=1; qGapBp+=L; q_ofs+=L
        elif op=='D':
            tGaps+=1; tGapBp+=L; t_ofs+=L
        # S/H 忽略或只影响 q_ofs=0
    q_aln=sum(blocks)
    mism = nm - (qGapBp + tGapBp)
    match = q_aln - mism
    return blocks, qS, tS, q_aln, q_aln, qGaps, qGapBp, tGaps, tGapBp, match, mism

# 3) 转换
with open(sam_path) as sf, open(psl_path,'w') as out:
    # 写 5 行空 header
    for _ in range(5):
        out.write('\t\n')

    for line in sf:
        if line.startswith('@'):
            continue
        flds = line.rstrip('\n').split('\t')
        qname = flds[0]
        flag  = int(flds[1])
        rname = flds[2]
        pos   = int(flds[3])
        cigar = flds[5]
        seq   = flds[9]
        # —— 新增：跳过不对齐记录 —— 
        if flag & 0x4:            # unmapped
            continue
        if cigar == '*' or cigar == '':
            continue

        # parse NM tag    
        nm = 0
        for tag in flds[11:]:
            if tag.startswith('NM:i:'):
                nm = int(tag.split(':')[-1])
                break

        blocks, qS, tS, q_aln, t_aln, qGaps, qGapBp, tGaps, tGapBp, match, mism = \
            parse_cigar(cigar, pos, nm)

        strand = '-' if (flag & 0x10) else '+'
        qsize  = len(seq)
        qstart = qS[0]
        qend   = qstart + q_aln
        tstart = pos - 1
        tend   = tstart + t_aln
        bc     = len(blocks)
        bs     = ','.join(map(str,blocks)) + ','
        qs     = ','.join(map(str,qS)) + ','
        ts     = ','.join(map(str,tS)) + ','

        row = [
            match, mism, 0, 0,
            qGaps, qGapBp, tGaps, tGapBp,
            strand, qname, qsize, qstart, qend,
            rname, chrom_size.get(rname,0), tstart, tend,
            bc, bs, qs, ts
        ]
        out.write('\t'.join(map(str,row)) + '\n')

print(f"Converted {sam_path} → {psl_path}")
