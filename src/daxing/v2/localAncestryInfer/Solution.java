package daxing.v2.localAncestryInfer;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import java.util.*;

public class Solution {

    SolutionElement[] solutionElements;
    int seqLen;
    Direction direction;

    public Solution(EnumSet<WindowSource.Source>[] candidateSolutionsSource, Direction direction){
        this.seqLen = candidateSolutionsSource.length;
        this.direction=direction;
        this.solutionElements=transform2Solution(candidateSolutionsSource, direction);
    }

    private Solution(SolutionElement[] solutionElements, int seqLen, Direction direction){
        this.solutionElements=solutionElements;
        this.seqLen=seqLen;
        this.direction=direction;
    }

    private SolutionElement[] transform2Solution(EnumSet<WindowSource.Source>[] candidateSolutions, Direction direction){
        if (direction==Direction.F){
            return transform2ForwardSolution(candidateSolutions);
        }else {
            return transform2ReverseSolution(candidateSolutions);
        }
    }

    private SolutionElement[] transform2ForwardSolution(EnumSet<WindowSource.Source>[] forwardCandidateSolutionsSource){
        EnumSet<WindowSource.Source> pastSources = forwardCandidateSolutionsSource[0];
        List<SolutionElement> solutionElementList = new ArrayList<>();
        SolutionElement currentSolutionElement=new SolutionElement(pastSources);
        int startIndex=0;
        int endIndex= 1;
        for (int i = 1; i < forwardCandidateSolutionsSource.length; i++) {
            if (pastSources.equals(forwardCandidateSolutionsSource[i])){
                endIndex++;
            }else {
                currentSolutionElement.setStartIndex(startIndex);
                currentSolutionElement.setEndIndex(endIndex);
                solutionElementList.add(currentSolutionElement);
                currentSolutionElement = new SolutionElement(forwardCandidateSolutionsSource[i]);
                pastSources=forwardCandidateSolutionsSource[i];
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

    private SolutionElement[] transform2ReverseSolution(EnumSet<WindowSource.Source>[] reverseCandidateSolutionsSource){
        EnumSet<WindowSource.Source> pastSources = reverseCandidateSolutionsSource[0];
        List<SolutionElement> solutionElementList = new ArrayList<>();
        SolutionElement currentSolutionElement=new SolutionElement(pastSources);
        int startIndex=0;
        int endIndex= 1;
        for (int i = 1; i < reverseCandidateSolutionsSource.length; i++) {
            if (pastSources.equals(reverseCandidateSolutionsSource[i])){
                endIndex++;
            }else {
                currentSolutionElement.setStartIndex(startIndex);
                currentSolutionElement.setEndIndex(endIndex);
                solutionElementList.add(currentSolutionElement);
                currentSolutionElement = new SolutionElement(reverseCandidateSolutionsSource[i]);
                pastSources=reverseCandidateSolutionsSource[i];
                startIndex = i;
                endIndex = i+1;
            }
        }
        currentSolutionElement.setStartIndex(startIndex);
        currentSolutionElement.setEndIndex(endIndex);
        solutionElementList.add(currentSolutionElement);

        SolutionElement[] solutionElements = new SolutionElement[solutionElementList.size()];
        EnumSet<WindowSource.Source> sources;
        SolutionElement solutionElement;
        int seqLen = solutionElementList.get(solutionElementList.size()-1).endIndex;
        for (int i = 0; i < solutionElements.length; i++) {
            sources = solutionElementList.get(i).getSources();
            solutionElement = new SolutionElement(sources);
            solutionElement.setStartIndex(seqLen-solutionElementList.get(i).endIndex);
            solutionElement.setEndIndex(seqLen-solutionElementList.get(i).startIndex);
            solutionElements[solutionElementList.size()-1-i] = solutionElement;
        }
        return solutionElements;
    }

    public SolutionElement[] getSolutionElements() {
        return solutionElements;
    }

    public Direction getDirection() {
        return direction;
    }

    public int getSeqLen() {
        return seqLen;
    }

    public int getSourceElementCount(){
        return this.solutionElements.length;
    }

    /**
     * 计算breakpoint的条件, ifSatisfiesCondition需为true
     * @return
     */
    public IntList getTargetSourceIndexList(){
        IntList targetSourceIndexList = new IntArrayList();
        EnumSet<WindowSource.Source> targetSources = EnumSet.range(WindowSource.Source.WE, WindowSource.Source.AT);
        WindowSource.Source source;
        for (int i = 0; i < this.getSourceElementCount(); i++) {
            if(this.solutionElements[i].getSources().size()==1){
                source =this.solutionElements[i].getSources().iterator().next();
                if (targetSources.contains(source)){
                    targetSourceIndexList.add(i);
                }
            }
        }
        return targetSourceIndexList;
    }

    /**
     * 交集
     * @param solutionElement target source对应的solution element
     */
    private void updateSolutionElement(int targetSourceIndex, SolutionElement solutionElement){
        if (this.solutionElements[targetSourceIndex].startIndex < solutionElement.startIndex){
            this.solutionElements[targetSourceIndex].startIndex=solutionElement.startIndex;
            if (targetSourceIndex > 0){
                this.solutionElements[targetSourceIndex-1].endIndex=solutionElement.startIndex;
            }
        }
        if (this.solutionElements[targetSourceIndex].endIndex > solutionElement.endIndex){
            this.solutionElements[targetSourceIndex].endIndex = solutionElement.endIndex;
            if (targetSourceIndex < this.getSourceElementCount()-1){
                this.solutionElements[targetSourceIndex+1].startIndex=solutionElement.endIndex;
            }
        }
    }

    public Solution extend(SolutionElement solutionElement){
        SolutionElement[] solutionElements = new SolutionElement[this.getSourceElementCount()+1];
        System.arraycopy(this.solutionElements, 0, solutionElements, 0, this.getSourceElementCount());
        solutionElements[this.getSourceElementCount()]=solutionElement;
        return new Solution(solutionElements, this.getSeqLen(), this.getDirection());
    }

    public static Solution coalescent(EnumSet<WindowSource.Source>[][] candidateSolutions, Direction direction){
        Solution solution, currentSolution;
        Set<Solution> solutionSet = new HashSet<>();
        List<Solution> optimumSolutionList = new ArrayList<>();
        EnumSet<WindowSource.Source> sources = EnumSet.range(WindowSource.Source.WE, WindowSource.Source.AT);
        for (int i = 0; i < candidateSolutions.length; i++) {
            solution = new Solution(candidateSolutions[i], direction);
            // Solution至少含有两种solution elements
            if (solution.getSourceElementCount()==1 && (!sources.contains(solution.getSolutionElements()[0].sources.iterator().next()))) continue;
            solutionSet.add(solution);
        }

        int miniCount=Integer.MAX_VALUE;
        for (Solution solution1: solutionSet){
            miniCount = solution1.getSourceElementCount() <  miniCount ? solution1.getSourceElementCount() : miniCount;
        }

        for (Solution elements: solutionSet){
            // 断点数目最少的solution
            if (elements.getSourceElementCount() != miniCount) continue;
            optimumSolutionList.add(elements);
        }

        solution = optimumSolutionList.get(0);
        IntList targetSourceIndexList0 = solution.getTargetSourceIndexList();
        IntList currentTargetSourceIndexList;
        for (int i = 1; i < optimumSolutionList.size(); i++) {
            currentSolution = optimumSolutionList.get(i);
            currentTargetSourceIndexList = currentSolution.getTargetSourceIndexList();
            if (!targetSourceIndexList0.equals(currentTargetSourceIndexList)) continue;
            for (int index : currentTargetSourceIndexList){
                solution.updateSolutionElement(index, currentSolution.getSolutionElements()[index]);
            }
        }
        return solution;
    }



    public static Solution calculateBreakPoint(EnumMap<Solution.Direction, EnumSet<WindowSource.Source>[][]> candidateSolutions){
        Solution forwardSolution = Solution.coalescent(candidateSolutions.get(Direction.F), Direction.F);
        EnumSet<WindowSource.Source> targetSourceEnum = EnumSet.range(WindowSource.Source.WE, WindowSource.Source.AT);
        int forwardCount, reverseCount;
        forwardCount = forwardSolution.getSourceElementCount();

        EnumSet<WindowSource.Source> forwardLastSourceEnum= forwardSolution.solutionElements[forwardCount-1].sources;
        if (forwardLastSourceEnum.size()==1 && targetSourceEnum.contains(forwardLastSourceEnum.iterator().next())){
            Solution reverseSolution = Solution.coalescent(candidateSolutions.get(Direction.R), Direction.R);
            reverseCount = reverseSolution.getSourceElementCount();
            if (forwardCount==1 || reverseCount==1) return forwardSolution;
            if (reverseSolution.getSolutionElements()[reverseCount-2].sources.equals(forwardLastSourceEnum)){
                if(!reverseSolution.getSolutionElements()[reverseCount-2].sources.equals(reverseSolution.getSolutionElements()[reverseCount-1])){
                    forwardSolution.updateSolutionElement(forwardCount-1, reverseSolution.getSolutionElements()[reverseCount-2]);
                    return forwardSolution.extend(reverseSolution.getSolutionElements()[reverseCount-1]);
                }
            }
        }
        return forwardSolution;
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

        public void setEndIndex(int endIndex){
            this.endIndex= endIndex;
        }

        public int getSourceNum(){
            return sources.size();
        }

        public int getPosCount(){
            return endIndex-startIndex;
        }

        @Override
        public boolean equals(Object o){
            if (o == this) return true;
            if (!(o instanceof SolutionElement)) return false;
            SolutionElement other = (SolutionElement)o;
            boolean sourcesEqual= (this.sources == null && other.sources == null) || (this.sources != null &&
                    this.sources.equals(other.sources));
            boolean startEqual = this.startIndex == other.startIndex;
            boolean endEqual = this.endIndex == other.endIndex;
            return sourcesEqual && startEqual && endEqual;
        }

        @Override
        public int hashCode() {
            return Objects.hash(sources, startIndex, endIndex);
        }
    }




}
