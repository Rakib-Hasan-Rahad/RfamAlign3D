-- MySQL dump 10.13  Distrib 8.0.41, for Linux (x86_64)
--
-- Host: mysql-rfam-live.ebi.ac.uk    Database: rfam_live
-- ------------------------------------------------------
-- Server version	5.7.27-30-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `motif_matches`
--

DROP TABLE IF EXISTS `motif_matches`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `motif_matches` (
  `motif_acc` varchar(7) NOT NULL,
  `rfam_acc` varchar(7) NOT NULL,
  `rfamseq_acc` varchar(25) DEFAULT NULL,
  `rfamseq_start` bigint(19) DEFAULT NULL,
  `rfamseq_stop` bigint(19) DEFAULT NULL,
  `query_start` int(11) DEFAULT NULL,
  `query_stop` int(11) DEFAULT NULL,
  `motif_start` int(11) DEFAULT NULL,
  `motif_stop` int(11) DEFAULT NULL,
  `e_value` varchar(15) DEFAULT NULL,
  `bit_score` double(7,2) DEFAULT NULL,
  KEY `motif_match_motif_acc_idx` (`motif_acc`),
  KEY `motif_match_rfam_acc_idx` (`rfam_acc`),
  KEY `motif_match_rfamseq_acc_idx` (`rfamseq_acc`),
  CONSTRAINT `motif_match_motif_acc` FOREIGN KEY (`motif_acc`) REFERENCES `motif_old` (`motif_acc`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `motif_match_rfam_acc` FOREIGN KEY (`rfam_acc`) REFERENCES `family` (`rfam_acc`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `motif_match_rfamseq_acc` FOREIGN KEY (`rfamseq_acc`) REFERENCES `rfamseq` (`rfamseq_acc`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1 ROW_FORMAT=COMPACT;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2026-01-07 17:33:47
