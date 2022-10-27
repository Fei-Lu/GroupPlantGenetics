package daxing.v2.localAncestryInfer;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import java.util.*;

public class Solution {

    // the first dim is forward and reverse, the second dim is SolutionElement
    EnumMap<Direction, SolutionElement[]> forwardReverseSolutionElementsMap;

    /**
     * Whether the start solution element satisfies the condition
     * condition: the start source do not belong to (WE || DE || FTT || AT)
     * 即开始source不能以WE, DE, FTT, AT开头
     */
    boolean ifSatisfiesCondition=false;
    int seqLen;

    public Solution(EnumSet<WindowSource.Source>[] forwardCandidateSolutionsSource,
                    EnumSet<WindowSource.Source>[] reverseCandidateSolutionsSource){
        this.seqLen = forwardCandidateSolutionsSource.length;
        this.forwardReverseSolutionElementsMap= new EnumMap<>(Direction.class);
        EnumSet<WindowSource.Source> fistForwardSolutionElement=forwardCandidateSolutionsSource[0];
        EnumSet<WindowSource.Source> fistReverseSolutionElement=forwardCandidateSolutionsSource[0];
        if (this.ifSatisfiesCondition(fistForwardSolutionElement) && this.ifSatisfiesCondition(fistReverseSolutionElement)){
            this.ifSatisfiesCondition=true;
        }
        SolutionElement[] solutionElements = this.transform2Solution(forwardCandidateSolutionsSource);
        this.forwardReverseSolutionElementsMap.put(Direction.F, solutionElements);
        solutionElements = this.transform2Solution(reverseCandidateSolutionsSource);
        this.forwardReverseSolutionElementsMap.put(Direction.R, solutionElements);
    }

    private SolutionElement[] transform2Solution(EnumSet<WindowSource.Source>[] candidateSolutionsSource){
        EnumSet<WindowSource.Source> pastSources = candidateSolutionsSource[0];
        List<SolutionElement> solutionElementList = new ArrayList<>();
        SolutionElement currentSolutionElement=new SolutionElement(pastSources);
        int startIndex=0;
        int endIndex= 1;
        for (int i = 1; i < candidateSolutionsSource.length; i++) {
            if (pastSources.equals(candidateSolutionsSource[i])){
                endIndex++;
            }else {
                currentSolutionElement.setStartIndex(startIndex);
                currentSolutionElement.setEndIndex(endIndex);
                solutionElementList.add(currentSolutionElement);
                currentSolutionElement = new SolutionElement(candidateSolutionsSource[i]);
                pastSources=candidateSolutionsSource[i];
                startIndex = i;
                endIndex = i+1;
            }
        }
        currentSolutionElement.setStartIndex(startIndex);
        currentSolutionElement.setEndIndex(endIndex);
        solutionElementList.add(currentSolutionElement);

        SolutionElement[] solutionElements = new SolutionElement[solutionElementList.size()];
        for (int i = 0; i < solutionElements.length; i++) {
            solutionElements[i] = solutionElementList.get(i);
        }
        return solutionElements;
    }

    private boolean ifSatisfiesCondition(EnumSet<WindowSource.Source> fistSolutionElement){
        if (fistSolutionElement.size()==1 && EnumSet.range(WindowSource.Source.WE, WindowSource.Source.AT).contains(fistSolutionElement.iterator().next())){
            return false;
        }
        return true;
    }

    public EnumMap<Direction, SolutionElement[]> getForwardReverseSolutionElementsMap() {
        return forwardReverseSolutionElementsMap;
    }

    public SolutionElement[] getSolutionElements(Direction direction){
        return this.getForwardReverseSolutionElementsMap().get(direction);
    }

    public int getSourceElementCount(Direction direction){
        return this.getForwardReverseSolutionElementsMap().get(direction).length;
    }

    /**
     * 计算breakpoint的条件, ifSatisfiesCondition需为true
     * @return
     */
    public IntList getTargetSourceIndexList(Direction direction){
        IntList targetSourceIndexList = new IntArrayList();
        EnumSet<WindowSource.Source> targetSources = EnumSet.range(WindowSource.Source.WE, WindowSource.Source.AT);
        WindowSource.Source targetSource;
        //  ifSatisfiesCondition == true, 即target source不会出现在i=0处
        for (int i = 1; i < this.getSourceElementCount(direction); i++) {
            if(this.getSolutionElements(direction)[i].getSourceNum()==1){
                targetSource = this.getSolutionElements(direction)[i].getSources().iterator().next();
                if (targetSources.contains(targetSource)){
                    targetSourceIndexList.add(i);
                }
            }
        }
        return targetSourceIndexList;
    }

    public enum Direction{
        F(0), R(1);

        int direction;

        Direction(int direction){
            this.direction=direction;
        }
    }

    public static class SolutionElement{
        EnumSet<WindowSource.Source> sources;
        int startIndex; // inclusive
        int endIndex; // exclusive

        public SolutionElement(EnumSet<WindowSource.Source> sources){
            this.sources= sources;
            this.startIndex = -1;
            this.endIndex = -1;
        }

        public EnumSet<WindowSource.Source> getSources() {
            return sources;
        }

        public void setStartIndex(int startIndex){
            this.startIndex = startIndex;
        }

        public void setEndIndex(int EndIndex){
            this.endIndex= endIndex;
        }

        public int getSourceNum(){
            return sources.size();
        }

        public int getPosCount(){
            return endIndex-startIndex;
        }

    }




}
