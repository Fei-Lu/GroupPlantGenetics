package daxing.v2.localAncestryInfer;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

public class Solution {

    SolutionElement[] solutionElements;

    public Solution(EnumSet<WindowSource.Source>[] candidateSolutionsSource){
        List<SolutionElement> solutionElementList = new ArrayList<>();
        EnumSet<WindowSource.Source> pastSources = candidateSolutionsSource[0];
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
        this.solutionElements=solutionElements;
    }

    public SolutionElement[] getSolutionElements() {
        return solutionElements;
    }

    public int getMaxFragmentCount(){
        return this.getSolutionElements().length;
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

        public void setStartIndex(int startIndex){
            this.startIndex = startIndex;
        }

        public void setEndIndex(int EndIndex){
            this.endIndex= endIndex;
        }
    }




}
