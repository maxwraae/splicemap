#!/usr/bin/python3

import sys,getopt,math

def bppt_get_pwm(pwmf):
    index = 0
    PWMBP = {}
    with open(pwmf,'r') as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith("#"):
                continue
            PWMBP[index] = {}
            line = tmp.split('\t')
            PWMBP[index]['A'] = float(line[1])
            PWMBP[index]['C'] = float(line[2])
            PWMBP[index]['G'] = float(line[3])
            PWMBP[index]['T'] = float(line[4])
            index += 1
    return(PWMBP)

def bppt_get_ppt(pptf):
    PPTS = {}
    with open(pptf,'r') as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith("#"):
                continue
            line = tmp.split('\t')
            PPTS[line[0]] = float(line[4])
    return(PPTS)

def bppt_bpscore(cbps):
    cbpsc = 1
    cbpsN = list(cbps)
    for i in range(0,len(cbpsN)):
        if i in PWMBP:
            cbpsc = cbpsc*PWMBP[i][cbpsN[i]]
        else:
            print("Error: the input bps is longer than the PWM")
    return(cbpsc)

def bppt_get_bpscore(l,basebp):
    nn = ['A','C','G','T']
    NNS = ['A','C','G','T']
    NN = []
    for i in range(2,l+1):
        for nns in NNS:
            for n in nn:
                newN = nns+n
                NN.append(newN)
        NNS = NN
        NN = []
    basebpsc = bppt_bpscore(basebp)
    cBPSC = {}
    for nn in NNS:
        cBPSC[nn] = bppt_bpscore(nn)/basebpsc
    return(cBPSC)

def bppt_get_AGEZ(seq,offset=12,maxL=-1):
    sL = len(seq)
    ss = seq[0:-offset].split("AG")
    if maxL > 0:
        pAG = sL - maxL
    else:
        pAG = sL - (offset+len(ss[-1])+14)
    if pAG < 0:
        pAG = 0
    return(pAG)

def bppt_get_pptsc(pptS,lppt,l_max_ppt,baseppt):
    end = len(pptS)-lppt
    if end > l_max_ppt:
        end = l_max_ppt
    pptbasesc = (l_max_ppt-lppt+1)*PPTS[baseppt]
    pptSC = 0
    for i in range(0,end):
        cppts = pptS[i:i+lppt]
        pptSC += PPTS[cppts]
    if pptbasesc == 0:
        return(pptSC)
    else:
        return((pptSC/pptbasesc))

def bppt_dis_pro(pAG,offset=22,interval=4):
    if abs(pAG-offset) > 700:
        return(0)
    else:
        return(1/(math.exp((abs(pAG-offset))/interval)+1))

def bppt_get_BPPTsc(seq,maxL,baseppt):
    lmotif = 7
    lppt = 8
    l_max_ppt = 20
    pstart = bppt_get_AGEZ(seq=seq,offset=12,maxL=maxL)
    sL = len(seq)
    totbpsc = 0
    totpptsc = 0
    totsc = 0
    npos = 0
    orinp = []
    for ipos in range(pstart,sL-3-lmotif):
        pAG = sL - ipos - 5
        bpS = seq[ipos:ipos+lmotif]
        bpSC = cBPSC[bpS]
        pptSC = 0
        dis3 = pAG - 1
        if dis3 > lppt + 3:
            pptS = seq[ipos+lmotif:sL-3]
            pptSC = bppt_get_pptsc(pptS,lppt,l_max_ppt,baseppt)
        SC = (bpSC * pptSC)
        inp = bpS+"\t"+str(pAG)+"\t"+str(bpSC)+"\t"+str(pptSC)+"\t"+str(SC)
        if len(orinp) == 0:
            orinp.append(inp)
        else:
            flag_in = 0
            for i in range(0,len(orinp)):
                line = orinp[i].split("\t")
                scold = float(line[4])
                if scold < SC:
                    orinp.insert(i,inp)
                    flag_in = 1
                    break
            if flag_in == 0:
                orinp.append(inp)
        totsc += SC
        totbpsc += bpSC
        totpptsc += pptSC
        npos += 1
    msc = totsc/npos
    mbpsc = totbpsc/npos
    mpptsc = totpptsc/npos
    dsc = []
    dbpsc = []
    dpptsc = []
    sdsc = 0
    sdbpsc = 0
    sdpptsc = 0
    for i in range(0,len(orinp)):
        line = orinp[i].split("\t")
        sc = float(line[4])
        bpsc = float(line[2])
        pptsc = float(line[3])
        dd = sc-msc
        dsc.append(dd)
        sdsc += dd*dd
        dd = bpsc-mbpsc
        dbpsc.append(dd)
        sdbpsc += dd*dd
        dd = pptsc-mpptsc
        dpptsc.append(dd)
        sdpptsc += dd*dd
    sdsc = math.sqrt(sdsc/npos)
    sdbpsc = math.sqrt(sdbpsc/npos)
    sdpptsc = math.sqrt(sdpptsc/npos)
    zsc = []
    zbps = []
    zppt = []
    for i in range(0,len(dsc)):
        zsc.append(dsc[i]/sdsc)
        zbps.append(dbpsc[i]/sdbpsc)
        zppt.append(dpptsc[i]/sdpptsc)
    return(orinp,zbps,zppt,zsc)

def bppt_print(idd,orinp,zbps,zppt,zsc,REPORTN):
    if REPORTN == 0:
        end = len(orinp)
    else:
        end = REPORTN
        if end > len(orinp):
            end = len(orinp)
    print("#id\tbps\tbp_pos\tsc_bps\tsc_ppt\tsc\tzsc_bps\tzsc_ppt\tzsc")
    for i in range(0,end):
        print(idd+"\t"+orinp[i]+"\t"+str(zbps[i])+"\t"+str(zppt[i])+"\t"+str(zsc[i]))

def main(argv):
    mL = 7
    basebp = "TACTAAC"
    baseppt = "TTTTTTTT"
    try:
        opts, args = getopt.getopt(argv,"b:r:p:i:h",["bppt"])
    except getopt.GetoptError as err:
        sys.exit(2)
    pwmf = pptf = fastaF = ""
    repn = 1
    for opt, arg in opts:
        if opt == '-h':
            print("Usage: BP_PPT.py -b -p -i -r -h")
            sys.exit()
        elif opt == '-b':
            pwmf = arg
        elif opt == '-p':
            pptf = arg
        elif opt == '-r':
            repn = arg
        elif opt == '-i':
            fastaF = arg
    REPORTN = int(repn)
    global PWMBP, PPTS, cBPSC
    PWMBP = bppt_get_pwm(pwmf)
    PPTS = bppt_get_ppt(pptf)
    cBPSC = bppt_get_bpscore(mL,basebp)
    idd = ""
    seq = ""
    with open(fastaF,'r') as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith(">"):
                if idd == "":
                    idd = tmp
                else:
                    (orinp,zbps,zppt,zsc) = bppt_get_BPPTsc(seq,maxL=-1,baseppt=baseppt)
                    bppt_print(idd,orinp,zbps,zppt,zsc,REPORTN)
                    idd = tmp
                    seq = ""
            else:
                seq = seq+tmp
    (orinp,zbps,zppt,zsc) = bppt_get_BPPTsc(seq,maxL=-1,baseppt=baseppt)
    bppt_print(idd,orinp,zbps,zppt,zsc,REPORTN)

if __name__ == "__main__":
    main(sys.argv[1:])
