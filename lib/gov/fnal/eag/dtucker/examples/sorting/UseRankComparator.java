package gov.fnal.eag.dtucker.examples.sorting;

//From http://java.sun.com/developer/JDCTechTips/2004/tt0716.html#1


import java.util.Collections;

public class UseRankComparator extends UseRank {

  void sort() {
    Collections.sort(miniDeck, new RankComparator());
  }

  public static void main(String[] args) {
    new UseRankComparator().exerciseDeck();
  }
}

