package abeel.genometools.gfa

import scala.collection.mutable.MutableList

case class FastSegment(val idx: Int, val sequenceLen: Int, var incomingCount: Int=0, var outgoingCount:Int=0) {
  override def toString() = {
    idx + "," + sequenceLen + "," + incomingCount + "," + outgoingCount
  }
}