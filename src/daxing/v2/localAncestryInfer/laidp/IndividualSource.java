package daxing.v2.localAncestryInfer.laidp;

import daxing.common.chrrange.ChrRange;
import daxing.common.utiles.IOTool;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

public class IndividualSource {

    String individualID;
    WindowSource[] windowSources;

    public IndividualSource(String fd_dxyFile, String individualID){
        this.windowSources=IndividualSource.getWindowSourceFrom(fd_dxyFile);
        this.individualID=individualID;
    }

    public static WindowSource[] getWindowSourceFrom(String fd_dxyFile){
        boolean containSimulate=fd_dxyFile.contains("simulate");
        int donorIndex = containSimulate ? 7 : 13;
        List<WindowSource> windowSourceList = new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(fd_dxyFile)) {
            String line;
            List<String> temp1, temp2;
            br.readLine();
            String refChr, currentChr=null;
            int refStart, refEnd, lastEnd = -1;
            String lastRefChr = null;
            EnumSet<WindowSource.Source> sourceEnumSet;
            WindowSource windowSource;
            ChrRange chrRange;
            temp1 = PStringUtils.fastSplit(br.readLine());
            boolean ifFirstLine=true;
            while ((line=br.readLine())!=null){
                temp2 = PStringUtils.fastSplit(line);
                if (ifFirstLine){
                    refChr = temp1.get(0);
                    currentChr = refChr;
                    refStart = Integer.parseInt(temp1.get(1));
                    refEnd = Integer.parseInt(temp2.get(1)) + 1;
                    chrRange = new ChrRange(refChr, refStart, refEnd);
                    sourceEnumSet = containSimulate ?
                            EnumSet.of(WindowSource.Source.getInstanceFromSubNum(Integer.parseInt(temp1.get(donorIndex))).get()) :
                            EnumSet.of(WindowSource.Source.valueOf(temp1.get(donorIndex)));
                    windowSource = new WindowSource(chrRange, sourceEnumSet);
                    windowSourceList.add(windowSource);
                    ifFirstLine=false;
                }
                refChr = temp2.get(0);
                if (refChr.equals(currentChr)){
                    refStart = Integer.parseInt(temp2.get(1));
                    refEnd = Integer.parseInt(temp1.get(2))+1;
                    chrRange = new ChrRange(refChr, refStart, refEnd);
                    sourceEnumSet = containSimulate ?
                            EnumSet.of(WindowSource.Source.getInstanceFromSubNum(Integer.parseInt(temp1.get(donorIndex))).get()
                                    , WindowSource.Source.getInstanceFromSubNum(Integer.parseInt(temp2.get(donorIndex))).get()):
                            EnumSet.of(WindowSource.Source.valueOf(temp1.get(donorIndex)),
                            WindowSource.Source.valueOf(temp2.get(donorIndex)));
                    windowSource = new WindowSource(chrRange, sourceEnumSet);
                    windowSourceList.add(windowSource);
                    lastRefChr = temp1.get(0);
                    lastEnd = Integer.parseInt(temp1.get(2));
                    temp1 = temp2;
                }else {
                    refStart = windowSourceList.get(windowSourceList.size()-1).getChrRange().getEnd();
                    refEnd = Integer.parseInt(temp1.get(2))+1;
                    chrRange = new ChrRange(currentChr, refStart, refEnd);
                    sourceEnumSet = containSimulate ?
                            EnumSet.of(WindowSource.Source.getInstanceFromSubNum(Integer.parseInt(temp1.get(donorIndex))).get()
                            , WindowSource.Source.getInstanceFromSubNum(Integer.parseInt(temp2.get(donorIndex))).get()):
                            EnumSet.of(WindowSource.Source.valueOf(temp1.get(donorIndex)),
                            WindowSource.Source.valueOf(temp2.get(donorIndex)));
                    windowSource = new WindowSource(chrRange, sourceEnumSet);
                    windowSourceList.add(windowSource);
                    temp1 = temp2;
                    ifFirstLine=true;
                    currentChr = refChr;
                }

            }
            chrRange = new ChrRange(lastRefChr, lastEnd+1, Integer.parseInt(temp1.get(2))+1);
            sourceEnumSet = containSimulate ?
                    EnumSet.of(WindowSource.Source.getInstanceFromSubNum(Integer.parseInt(temp1.get(donorIndex))).get()):
                    EnumSet.of(WindowSource.Source.valueOf(temp1.get(donorIndex)));
            windowSource = new WindowSource(chrRange, sourceEnumSet);
            windowSourceList.add(windowSource);
        } catch (IOException e) {
            e.printStackTrace();
        }
        Collections.sort(windowSourceList);
        return windowSourceList.toArray(new WindowSource[windowSourceList.size()]);
    }

    public WindowSource[] getWindowSources() {
        return windowSources;
    }

    public String[] getChr(){
        Set<String> chrSet = new HashSet<>();
        for (int i = 0; i < windowSources.length; i++) {
            chrSet.add(windowSources[i].getChrRange().getChr());
        }
        List<String> chrsList = new ArrayList<>(chrSet);
        Collections.sort(chrsList);
        String[] res = new String[chrsList.size()];
        for (int i = 0; i < chrsList.size(); i++) {
            res[i]=chrsList.get(i);
        }
        return res;
    }

    public WindowSource[] getSubsetWindowSources(String refChr){
        WindowSource query =  new WindowSource(new ChrRange(refChr, -2, -1), null);
        int chrStartHit = Arrays.binarySearch(this.getWindowSources(), query);
        int chrStartIndex = chrStartHit < 0 ? -chrStartHit-1: chrStartHit;
        query = new WindowSource(new ChrRange(refChr, Integer.MAX_VALUE, Integer.MAX_VALUE), null);
        int chrEndHit = Arrays.binarySearch(this.getWindowSources(), query);
        int chrEndIndex = chrEndHit < 0 ? -chrEndHit-2: chrEndHit;
        WindowSource[] res = new WindowSource[chrEndIndex-chrStartIndex+1];
        System.arraycopy(this.getWindowSources(), chrStartIndex, res, 0, chrEndIndex-chrStartIndex+1);
        return res;
    }

    public String getIndividualID() {
        return individualID;
    }

    /**
     * return index of windows which have inconsistent ancestry
     */
    public TIntArrayList getInconsistentAncestry(){
        TIntArrayList indexList=new TIntArrayList();
        for (int i = 0; i < windowSources.length; i++) {
            if (windowSources[i].sources.size() > 1){
                indexList.add(i);
            }
        }
        return indexList;
    }

    /**
     * Determination of local continuous introgression windows
     * @param conjunctionNum local
     * @return
     */
    public WindowSource[] selectCandidateWindow(int conjunctionNum){
        List<WindowSource> windowSourceList = new ArrayList<>();
        String[] chrs = this.getChr();
        WindowSource windowSource;
        for (int j = 0; j < chrs.length; j++) {
            // introgression window index set
            WindowSource[] windowSources = this.getSubsetWindowSources(chrs[j]);
            TIntSet introgressionWindowIndexSet = new TIntHashSet();
            for (int i = 0; i < windowSources.length; i++) {
                if (windowSources[i].containAnySourceIn(EnumSet.range(WindowSource.Source.WE, WindowSource.Source.AT))){
                    introgressionWindowIndexSet.add(i);
                }
            }
            TIntSet resultIndexSet =  new TIntHashSet();
            TIntIterator intIterator = introgressionWindowIndexSet.iterator();
            while (intIterator.hasNext()){
                int index = intIterator.next();
                if (index > 1 && index < windowSources.length-conjunctionNum){
                    int i=1;
                    while (i <= conjunctionNum){
                        resultIndexSet.add(index-i < 0 ? 0:index-i);
                        resultIndexSet.add(index+i);
                        i++;
                    }

                }else if (index == 1){
                    resultIndexSet.add(index-1);
                }else if (index == (windowSources.length-conjunctionNum)){
                    int i = 1;
                    while (i < conjunctionNum){
                        resultIndexSet.add(index+i);
                        i++;
                    }
                }else if (index == 0){
                    int i = 0;
                    while (i < conjunctionNum){
                        resultIndexSet.add(i);
                        i++;
                    }
                }
            }

            if (resultIndexSet.size()==0) continue;

            // introgression window index set to sorted arrat
            int[] resultIndexArray = resultIndexSet.toArray();
            Arrays.sort(resultIndexArray);

            // int[2] start end represent successive introgression window index
            List<int[]> indexList = new ArrayList<>();
            int initializeStartIndex=resultIndexArray[0];
            int initializeEndIndex=initializeStartIndex;
            int[] startEnd;
            for (int i = 1; i < resultIndexArray.length; i++) {
                if (resultIndexArray[i]-resultIndexArray[i-1]==1){
                    initializeEndIndex = resultIndexArray[i];
                }else {
                    startEnd = new int[2];
                    startEnd[0] = initializeStartIndex;
                    startEnd[1] = initializeEndIndex;
                    indexList.add(startEnd);
                    initializeStartIndex = resultIndexArray[i];
                    initializeEndIndex = initializeStartIndex ;
                }
            }
            startEnd = new int[2];
            startEnd[0]=initializeStartIndex;
            startEnd[1] = initializeEndIndex;
            indexList.add(startEnd);

            // add ChrRange to chrRangeList
            ChrRange chrRange;
            int refStart, refEnd;
            Set<WindowSource.Source> sourceSet;
            EnumSet<WindowSource.Source> sourceEnumSet;
            for (int i = 0; i < indexList.size(); i++) {
                startEnd = indexList.get(i);
                if (startEnd[0] == startEnd[1]) continue;
                sourceSet = new HashSet<>();
                for (int k = startEnd[0]; k < startEnd[1]; k++) {
                    sourceSet.addAll(windowSources[k].getSourceSet());
                }
                refStart=windowSources[startEnd[0]].getChrRange().getStart();
                refEnd = windowSources[startEnd[1]].getChrRange().getEnd();
                chrRange = new ChrRange(chrs[j], refStart, refEnd);
                sourceEnumSet=EnumSet.copyOf(sourceSet);
                windowSource = new WindowSource(chrRange, sourceEnumSet);
                windowSourceList.add(windowSource);
            }
        }
        WindowSource[] windowSources = new WindowSource[windowSourceList.size()];
        for (int i = 0; i < windowSourceList.size(); i++) {
            windowSources[i]=windowSourceList.get(i);
        }
        Arrays.sort(windowSources);
        return windowSources;
    }

    private int[] getStartPos(String refChr){
        WindowSource[] refChrWindowSource = this.getSubsetWindowSources(refChr);
        int[] startPoss=new int[refChrWindowSource.length];
        for (int i = 0; i < startPoss.length; i++) {
            startPoss[i]=refChrWindowSource[i].getChrRange().getStart();
        }
        return startPoss;
    }

    private int[] getEndPos(String refChr){
        WindowSource[] refChrWindowSource = this.getSubsetWindowSources(refChr);
        int[] endPoss = new int[refChrWindowSource.length];
        for (int i = 0; i < endPoss.length; i++) {
            endPoss[i]= refChrWindowSource[i].getChrRange().getEnd();
        }
        return endPoss;
    }

    /**
     * start target source指的是在一个windowSource区间内首次出现的非NONE的source
     * @param windowSource
     * @return start target source
     */
    public EnumSet<WindowSource.Source> getStartTargetSourceFrom(WindowSource windowSource){
        String refChr = windowSource.getChrRange().getChr();
        int[] startPoss = this.getStartPos(refChr);
        int[] endPoss = this.getEndPos(refChr);
        int queryStart=windowSource.getChrRange().getStart();
        int queryEnd = windowSource.getChrRange().getEnd();
        int startIndex = Arrays.binarySearch(startPoss, queryStart);
        int endIndex = Arrays.binarySearch(endPoss, queryEnd);
        WindowSource[] refChrWindowSource = this.getSubsetWindowSources(refChr);
        EnumSet<WindowSource.Source> sources = EnumSet.range(WindowSource.Source.WE, WindowSource.Source.AT);
        for (int i = startIndex; i <= endIndex; i++) {
            if (refChrWindowSource[i].containAnySourceIn(sources)){
                return refChrWindowSource[i].getSources();
            }
        }
        return null;
    }

    /**
     * end target source指的是在一个windowSource区间内最后一次出现的非NONE的source
     * @param windowSource
     * @return end target source
     */
    public EnumSet<WindowSource.Source> getEndTargetSourceFrom(WindowSource windowSource){
        String refChr = windowSource.getChrRange().getChr();
        int[] startPoss = this.getStartPos(refChr);
        int[] endPoss = this.getEndPos(refChr);
        int queryStart=windowSource.getChrRange().getStart();
        int queryEnd = windowSource.getChrRange().getEnd();
        int startIndex = Arrays.binarySearch(startPoss, queryStart);
        int endIndex = Arrays.binarySearch(endPoss, queryEnd);
        WindowSource[] refChrWindowSource = this.getSubsetWindowSources(refChr);
        EnumSet<WindowSource.Source> sources = EnumSet.range(WindowSource.Source.WE, WindowSource.Source.AT);
        for (int i = endIndex; i >= startIndex; i--) {
            if (refChrWindowSource[i].containAnySourceIn(sources)){
                return refChrWindowSource[i].getSources();
            }
        }
        return null;
    }


    public static void test(String genotypeFile, String fd_dxyFile, String groupByPop2IndividualFile, String outFile){
        String individualID="LR_0002";
//        GenoGrid genoGrid = new GenoGrid(genotypeFile, GenoIOFormat.VCF_GZ);
        IndividualSource individualSource = new IndividualSource(fd_dxyFile, individualID);
        TIntArrayList indexList = individualSource.getInconsistentAncestry();
        System.out.println();
    }

}
