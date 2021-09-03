package strsCluster;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.io.FileReader;
import java.util.regex.*;
import java.util.List;

public class exCDhit {

    /**
     * 提取cdhit聚类文件的聚类信息，每个类一个list
     * @param path
     * @return names
     * @throws IOException
     */
    public String[][] readClstr (String path) throws IOException {
        List<String[]> clstrs = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String temp, target = null;
            List<String> names = new ArrayList<>();
            Pattern p = Pattern.compile("(\\d+).*(>\\w+).*\\s+(\\*?).*");
            while ((temp = br.readLine()) != null) {
                if (temp.charAt(0) == '>') {
                    if (target != null) {
                        String[] cone = new String[names.size() + 1];
                        cone[0] = target.toLowerCase();
                        for (int i = 1; i < cone.length; i++) cone[i] = names.get(i-1).toLowerCase();
                        clstrs.add(cone);
                        target = null; names = new ArrayList<>();
                    }
                }
                else {
                    Matcher m = p.matcher(temp);
                    if (m.matches()) {
                        if (m.group(3).equals("*")) target = m.group(2);
                        else names.add(m.group(2));
                    }   
                }
            }
            if (target != null) {
                String[] cone = new String[names.size() + 1];
                cone[0] = target.toLowerCase();
                for (int i = 1; i < cone.length; i++) cone[i] = names.get(i-1).toLowerCase();
                clstrs.add(cone);
            }
        }
        return clstrs.toArray(String[][]::new);
    }
}
