function clarified_string = conditionStringToFilename(condition_string)
    proxy_string = condition_string;
    proxy_string = strrep(proxy_string, 'STANCE_LEFT', 'stanceL');
    proxy_string = strrep(proxy_string, 'STANCE_RIGHT', 'stanceR');
    proxy_string = strrep(proxy_string, 'ILLUSION_LEFT', 'illuL');
    proxy_string = strrep(proxy_string, 'ILLUSION_RIGHT', 'illuR');
    proxy_string = strrep(proxy_string, 'ONE', 'step1');
    proxy_string = strrep(proxy_string, 'TWO', 'step2');
    proxy_string = strrep(proxy_string, 'THREE', 'step3');
    proxy_string = strrep(proxy_string, 'FOUR', 'step4');
    clarified_string = proxy_string;
end